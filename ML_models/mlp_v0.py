import pandas as pd
import torch
from torch import nn
from torch.utils.data import TensorDataset, random_split, DataLoader

import pytorch_lightning as pl
from pytorch_lightning import LightningModule
from pytorch_lightning.loggers import TensorBoardLogger
from pytorch_lightning.callbacks.early_stopping import EarlyStopping

from matplotlib import pyplot as plt
import seaborn as sns
sns.set()

seeds = [159765, 1597652534, 784642, 19807, 10, 389471, 182, 131, 907078]
class MLP(LightningModule):
    def __init__(self, args, seed=seeds[5]):
        super().__init__()
        self.hparams.seed = seed
        pl.seed_everything(self.hparams.seed)
        self.save_hyperparameters(args)

        self.loss_fn = nn.L1Loss()
        self.act_fn = nn.ReLU()

        n = self.hparams.nodes_hlayers
        d = self.hparams.drop_vals

        self.layers = nn.Sequential(

            nn.Linear(self.hparams.input_size, n[0]),
            self.act_fn,
            #nn.BatchNorm1d(n[0]),
            # nn.Dropout(p=d[0]),

            nn.Linear(n[0], n[1]),
            self.act_fn,
           # nn.BatchNorm1d(n[1]),
            # nn.Dropout(p=d[1]),

            nn.Linear(n[1], n[2]),
            self.act_fn,
         #   nn.BatchNorm1d(n[2]),
            # nn.Dropout(p=d[2]),

            nn.Linear(n[2], n[3]),
            self.act_fn,
        #    nn.BatchNorm1d(n[3]),
            # nn.Dropout(p=d[3]),

            nn.Linear(n[3], 1),
        )

        self.kaiming_normal('relu')

        self.y = []
        self.ypred = []

    def kaiming_normal(self, nl):
        for param in self.parameters():
            if len(param.shape) == 1:
                nn.init.constant_(param, 0.1)
            else:
                nn.init.kaiming_normal_(param, mode='fan_in', nonlinearity=nl)
        self.hparams.kaiming_normal = 1

    def forward(self, x):
        return self.layers(x)

    def configure_optimizers(self):
        return torch.optim.SGD(self.parameters(), lr=self.hparams.learning_rate, momentum=0.9, nesterov=1)

    def training_step(self, train_batch, batch_idx):
        x, y = train_batch
        y_pred = torch.flatten(self.forward(x))
        #y_pred, y = self.unscale('s_bias', y_pred), self.unscale('s_bias', y)
        loss = self.loss_fn(y_pred, y)

        self.log('train_loss', loss)

        return {'loss': loss}

    def training_epoch_end(self, outputs):
        avg_loss = torch.stack([x['loss'] for x in outputs]).mean()
        self.logger.experiment.add_scalars('avg_loss', {'train': avg_loss}, self.current_epoch)

    def validation_step(self, valid_batch, batch_idx):
        x, y = valid_batch
        y_pred = torch.flatten(self.forward(x))
        #y_pred, y = self.unscale('s_bias', y_pred), self.unscale('s_bias', y)
        loss = self.loss_fn(y_pred, y)

        self.log('val_loss', loss, prog_bar=1)

        return {'loss': loss}

    def validation_epoch_end(self, outputs):
        avg_loss = torch.stack([x['loss'] for x in outputs]).mean()
        self.logger.experiment.add_scalars('avg_loss', {'val': avg_loss}, self.current_epoch)

    def test_step(self, batch, batch_idx):
        x, y = batch
        y_pred = torch.flatten(self.forward(x))
        #y_pred, y = self.unscale('s_bias', y_pred), self.unscale('s_bias', y)
        loss = self.loss_fn(y_pred, y)

        self.log('test_loss', loss)

        self.y = y.tolist()
        self.ypred = y_pred.tolist()

    def unscale(self, name, data):
        mean, std = self.hparams.distinfo[name]

        return data * std + mean

    def setup(self, stage=None):

        if stage == 'fit':
            df = pd.read_csv(self.hparams.csvdata_path).drop('nbonds', axis=1)

            df = df.iloc[:, :df.columns.get_loc('s_bias') + 1]
            feature_names = df.columns[:-1]
            self.standardize(df, stage)

            target = torch.tensor(df['s_bias'].values).float()

            features = torch.tensor(df.drop('s_bias', axis=1).values).float()
            tensordata = TensorDataset(features, target)
            data_size = len(tensordata)

            num_val = int(self.hparams.val_frac * data_size)
            num_train = data_size - num_val
            self.train_set, self.val_set = random_split(tensordata, (num_train, num_val))

        if stage == 'test':
            df = pd.read_csv('crambin_data.csv').drop('nbonds', axis=1)

            df = df.iloc[:, :df.columns.get_loc('s_bias') + 1]
            self.standardize(df, stage)

            target = torch.tensor(df['s_bias'].values).float()
            features = torch.tensor(df.drop('s_bias', axis=1).values).float()
            self.test_set = TensorDataset(features, target)

    def standardize(self, df, stage):
        if stage != 'test':
            self.hparams.distinfo = dict()
            for el in df.columns[:-1]:
                mean, std = df[el].mean(), df[el].std()

                if std > 0:
                    df[el] = (df[el] - mean) / std

                self.hparams.distinfo[el] = (mean, std)

        else:
            for el in df.columns[:-1]:
                mean, std = self.hparams.distinfo[el][0], self.hparams.distinfo[el][1]

                if std > 0:
                    df[el] = (df[el] - mean) / std

    def train_dataloader(self):
        return DataLoader(self.train_set, batch_size=self.hparams.batch_size, shuffle=True)

    def val_dataloader(self):
        return DataLoader(self.val_set, batch_size=self.hparams.batch_size)

    def test_dataloader(self):
        return DataLoader(self.test_set, batch_size=1023)


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument('-lr', '--learning_rate', type=float, default=0.000001)
    parser.add_argument('-lrf', '--lr_find', action="store_true", default=0)
    parser.add_argument('-nls', '--nodes_hlayers', type=tuple, default=(128, 64, 32, 16))
    parser.add_argument('-dv', '--drop_vals', type=tuple, default=(0.0, 0.0, 0.0, 0.0))

    parser.add_argument('-csv', '--csvdata_path', type=str, default='./mlp_data_v6.csv')
    parser.add_argument('-isz', '--input_size', type=int, default=7)
    parser.add_argument('-bsz', '--batch_size', type=int, default=64)
    parser.add_argument('-vf', '--val_frac', type=float, default=0.2)

    parser.add_argument('--gpus', type=int, default=0)

    args = parser.parse_args()
    model = MLP(args)

    test = 1
    version_num = 12
    name = 'test_logs' if test else f'mlp_logs_v{version_num}'

    logger = TensorBoardLogger('light_logs', name=name)
    early_stop_callback = EarlyStopping(monitor="val_loss", min_delta=0.002, patience=10, verbose=1, mode="min")
    trainer = pl.Trainer(max_epochs=250, gpus=args.gpus, logger=logger,
                         callbacks=[early_stop_callback], stochastic_weight_avg=1)

    if args.lr_find:
        lr_finder = trainer.tuner.lr_find(model)
        lr_finder.plot(suggest=1, show=1)
        args.learning_rate = lr_finder.suggestion()
        model = MLP(args)

    trainer.fit(model)
    trainer.test(model)

    sns.set()
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, dpi=300)
    fig.set_size_inches(10, 10)

    inds = range(len(model.y))
    residuals = [model.ypred[i] - model.y[i] for i in inds]

    sns.scatterplot(x=inds, y=model.y, ax=ax1, label='Actual')
    sns.scatterplot(x=inds, y=model.ypred, ax=ax2, label='Prediction')
    sns.scatterplot(x=inds, y=model.y, ax=ax3, label='Actual')
    sns.scatterplot(x=inds, y=model.ypred, ax=ax3, label='Prediction')

    for ax in (ax1, ax2, ax3):
        ax.set(xlabel='row_num', ylabel='s_bias')
        ax.legend()

    fig.show()
    fig.savefig(f'{logger.save_dir}/{logger.name}/version_{logger.version}/sbias_preds.png')

    fig, (ax1, ax2) = plt.subplots(2, 1, dpi=300)
    fig.set_size_inches(6, 8)

    sns.scatterplot(x=inds, y=residuals, ax=ax1)
    sns.histplot(data=residuals, ax=ax2)

    ax1.set(xlabel='Row number', ylabel='Residual', title='Residual Plot')
    ax2.set(xlabel='Residual', ylabel='Count', title='Residual Distribution')

    fig.show()
    fig.savefig(f'{logger.save_dir}/{logger.name}/version_{logger.version}/residual_info.png')