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

def RMSELoss(ypred, y):
    return ((ypred - y) ** 2).mean().sqrt()

class MLP(LightningModule):
    def __init__(self, args):
        super().__init__()
        self.save_hyperparameters(args)
        pl.seed_everything(self.hparams.seed)

        acts = {'relu': nn.ReLU(), 'leaky_relu': nn.LeakyReLU(), 'selu': nn.SELU()}
        loss_fns = {'L1': nn.L1Loss(), 'smooth_L1': nn.SmoothL1Loss(), 'MSE': nn.MSELoss(), 'RMSE': RMSELoss}

        self.act_fn = acts[self.hparams.act]
        self.loss_fn = loss_fns[self.hparams.loss]

        n = self.hparams.nodes_hlayers

        self.lin_layers = nn.ModuleList(
            [nn.Linear(self.hparams.input_size, n[0])] +
            [nn.Linear(n[i], n[i + 1]) for i in range(len(n) - 1)] +
            [nn.Linear(n[-1], 1)])

        if self.hparams.batchnorm:
            self.bn_layers = nn.ModuleList([nn.BatchNorm1d(el) for el in n])

        if self.hparams.dropout:
            self.do_layers = nn.ModuleList([nn.Dropout(el) for el in self.hparams.drop_vals])

        if self.hparams.kaiming_normal[0]:
            self.kaiming_normal(self.hparams.kaiming_normal[1])

        self.y = []
        self.ypred = []

    def kaiming_normal(self, nl):
        for param in self.parameters():
            if len(param.shape) == 1:
                nn.init.constant_(param, 0.1)
            else:
                nn.init.kaiming_normal_(param, mode='fan_in', nonlinearity=nl)

    def forward(self, x):
        for i in range(len(self.lin_layers) - 1):
            x = self.act_fn(self.lin_layers[i](x))

            if self.hparams.batchnorm:
                x = self.bn_layers[i](x)

            if self.hparams.dropout:
                x = self.do_layers[i](x)

        return self.lin_layers[-1](x)

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
            df = pd.read_csv(self.hparams.csvdata_path).drop(['structure_id'], axis=1)
            df = df.iloc[:, :df.columns.get_loc('s_bias') + 1]

            #self.standardize(df, stage)

            target = torch.tensor(df['s_bias'].values).float()

            features = torch.tensor(df.drop('s_bias', axis=1).values).float()
            tensordata = TensorDataset(features, target)
            data_size = len(tensordata)

            num_val = int(self.hparams.val_frac * data_size)
            num_train = data_size - num_val
            self.train_set, self.val_set = random_split(tensordata, (num_train, num_val))

        if stage == 'test':
            df = pd.read_csv('crambin_data.csv').drop(['structure_id'], axis=1)

            df = df.iloc[:, :df.columns.get_loc('s_bias') + 1]
           # self.standardize(df, stage)

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

    parser.add_argument('-lr', '--learning_rate', type=float, default=1e-5)
    parser.add_argument('--seed', type=int, default=389470)
    parser.add_argument('-lrf', '--lr_find', action="store_true", default=0)
    parser.add_argument('-nls', '--nodes_hlayers', type=tuple, default=(64, 32, 16, 8))
    parser.add_argument('-dv', '--drop_vals', type=tuple, default=(0.2, 0.2, 0.2, 0.0))
    parser.add_argument('-bn', '--batchnorm', action="store_true", default=1)
    parser.add_argument('-do', '--dropout', action="store_true", default=0)

    parser.add_argument('-af', '--act', type=str, default='leaky_relu')
    parser.add_argument('-lf', '--loss', type=str, default='MSE')
    parser.add_argument('-kn', '--kaiming_normal', type=tuple, default=(1, 'leaky_relu'))

    parser.add_argument('-csv', '--csvdata_path', type=str, default='./mlp_data_v6.csv')
    parser.add_argument('-isz', '--input_size', type=int, default=8)
    parser.add_argument('-bsz', '--batch_size', type=int, default=32)
    parser.add_argument('-vf', '--val_frac', type=float, default=0.2)

    parser.add_argument('--gpus', type=int, default=0)

    args = parser.parse_args()
    model = MLP(args)

    test = 1
    version_num = 12
    name = 'test_logs' if test else f'mlp_logs_v{version_num}'

    logger = TensorBoardLogger('light_logs', name=name)
    early_stop_callback = EarlyStopping(monitor="val_loss", min_delta=0.08, patience=8, verbose=1, mode="min")
    trainer = pl.Trainer(max_epochs=500, gpus=args.gpus, logger=logger,
                         callbacks=[early_stop_callback], stochastic_weight_avg=1)

    if args.lr_find:
        lr_finder = trainer.tuner.lr_find(model)
        lr_finder.plot(suggest=1, show=1)
        args.learning_rate = lr_finder.suggestion()
        model = MLP(args)

    trainer.fit(model)
    trainer.test(model)

    sns.set(style='white')

    # Plot actual and predicted entropy values
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, dpi=300)
    fig.set_size_inches(10, 10)

    inds = range(len(model.y))

    sns.scatterplot(x=inds, y=model.y, ax=ax1, label='Actual')
    sns.scatterplot(x=inds, y=model.ypred, ax=ax2, label='Prediction')
    sns.scatterplot(x=inds, y=model.y, ax=ax3, label='Actual')
    sns.scatterplot(x=inds, y=model.ypred, ax=ax3, label='Prediction')

    for ax in (ax1, ax2, ax3):
        ax.legend()

    fig.suptitle('Actual & Predicted Configurational Entropy')
    fig.supxlabel('Row Number')
    fig.supylabel('Entropy')
    fig.show()
    fig.savefig(f'{logger.save_dir}/{logger.name}/version_{logger.version}/sbias_preds.png')


    # Find and plot the residuals for the fit
    residuals = [model.ypred[i] - model.y[i] for i in inds]
    fig, (ax1, ax2) = plt.subplots(2, 1, dpi=300)
    fig.set_size_inches(6, 8)

    sns.scatterplot(x=inds, y=residuals, ax=ax1)
    sns.histplot(data=residuals, ax=ax2)

    ax1.set(xlabel='Row number', ylabel='Residual', title='Residual Plot')
    ax2.set(xlabel='Residual', ylabel='Count', title='Residual Distribution')

    fig.show()
    fig.savefig(f'{logger.save_dir}/{logger.name}/version_{logger.version}/residual_info.png')


    # Find and plot the percentage errors for the fit
    pe = [100 * (model.ypred[i] - model.y[i]) / model.y[i] for i in inds]
    fig, (ax1, ax2) = plt.subplots(2, 1, dpi=300)
    fig.set_size_inches(6, 8)

    sns.scatterplot(x=inds, y=pe, ax=ax1)
    sns.histplot(data=pe, ax=ax2)

    ax1.set(xlabel='Row number', ylabel='Percent error', title='Percent Error Plot')
    ax2.set(xlabel='Percent error', ylabel='Count', title='Percent Error Distribution')

    fig.show()
    fig.savefig(f'{logger.save_dir}/{logger.name}/version_{logger.version}/pe_info.png')



