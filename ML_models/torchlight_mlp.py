import pandas as pd
import torch
from torch import nn
from torch.utils.data import TensorDataset, random_split, DataLoader

import pytorch_lightning as pl
from pytorch_lightning import LightningModule, LightningDataModule
from pytorch_lightning.loggers import TensorBoardLogger
from pytorch_lightning.callbacks.early_stopping import EarlyStopping


class MLP(LightningModule):
    def __init__(self, data_hparams, learn_rate=0.01, numnodes_hiddenlayers=(12, 10, 7), dropout_vals=(0.5, 0.4, 0.2)):
        super().__init__()

        self.dhp = data_hparams
        self.lr = learn_rate
        self.n_hl = numnodes_hiddenlayers
        self.dropv = dropout_vals
        self.loss_fn = nn.MSELoss()

        n1, n2, n3 = self.n_hl
        d1, d2, d3 = self.dropv

        self.layers = nn.Sequential(

            nn.Linear(self.dhp['input_size'], n1),
            nn.BatchNorm1d(n1),
            nn.Dropout(d1),
            nn.ReLU(),

            nn.Linear(n1, n2),
            nn.BatchNorm1d(n2),
            nn.Dropout(d2),
            nn.ReLU(),

            nn.Linear(n2, n3),
            nn.BatchNorm1d(n3),
            nn.Dropout(d3),
            nn.ReLU(),

            nn.Linear(n3, 1),
        )

        self.save_hyperparameters()

    def forward(self, x):
        return self.layers(x)

    def configure_optimizers(self):
        return torch.optim.Adam(self.parameters(), lr=self.lr)

    def training_step(self, train_batch, batch_idx):
        x, y = train_batch
        y_pred = torch.flatten(self.forward(x))
        loss = self.loss_fn(y_pred, y)

        self.log('train_loss', loss)
        self.logger.experiment.add_scalars('loss', {'train': loss}, self.global_step)

        return {'loss': loss}

    def training_epoch_end(self, outputs):
        avg_loss = torch.stack([x['loss'] for x in outputs]).mean()
        self.logger.experiment.add_scalars('avg_loss', {'train': avg_loss}, self.current_epoch)

    def validation_step(self, valid_batch, batch_idx):
        x, y = valid_batch
        y_pred = torch.flatten(self.forward(x))
        loss = self.loss_fn(y_pred, y)

        self.log('val_loss', loss, prog_bar=1)
        self.logger.experiment.add_scalars('loss', {'val': loss}, self.global_step)

        return {'loss': loss}

    def validation_epoch_end(self, outputs):
        avg_loss = torch.stack([x['loss'] for x in outputs]).mean()
        self.logger.experiment.add_scalars('avg_loss', {'val': avg_loss}, self.current_epoch)

    def test_step(self, batch, batch_idx):
        x, y = batch
        y_pred = torch.flatten(self.forward(x))
        loss_fn = nn.MSELoss()
        loss = loss_fn(y_pred, y)

        self.log('test_loss', loss)


class MLPDataMod(LightningDataModule):

    def __init__(self, data_dir='./mlp_data.csv', batch_size=1000, val_frac=0.2):
        super().__init__()
        self.data_dir = data_dir
        self.batch_size = batch_size
        self.val_frac = val_frac

    def setup(self, stage=None):
        self.df = pd.read_csv(self.data_dir)

        feature_names = self.df.columns[:-1]
        self.input_size = len(feature_names)
        self.standardize(feature_names)

        target = torch.tensor(self.df['s_bias'].values).float()
        features = torch.tensor(self.df.drop('s_bias', axis=1).values).float()
        tensordata = TensorDataset(features, target)
        data_size = len(tensordata)

        num_val = int(self.val_frac * data_size)
        num_train = data_size - num_val

        self.train_set, self.val_set = random_split(tensordata, (num_train, num_val))

    def standardize(self, feature_names):
        for el in feature_names:
            mean, std = self.df[el].mean(), self.df[el].std()
            self.df[el] = (self.df[el] - mean) / std

    def train_dataloader(self):
        return DataLoader(self.train_set, batch_size=self.batch_size, shuffle=True)

    def val_dataloader(self):
        return DataLoader(self.val_set, batch_size=self.batch_size)


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument('-lr', '--learn_rate', type=float, default=0.01)
    parser.add_argument('-lrf', '--lr_find', action="store_true", default=False)
    parser.add_argument('-nls', '--nodes_hlayers', type=tuple, default=(12, 10, 7))
    parser.add_argument('-dv', '--drop_vals', type=tuple, default=(0.5, 0.4, 0.2))

    parser.add_argument('-dpt', '--data_path', type=str, default='./mlp_data_v2.csv')
    parser.add_argument('-bsz', '--batch_size', type=int, default=10)
    parser.add_argument('-vf', '--val_frac', type=float, default=0.2)

    args = parser.parse_args()

    logger = TensorBoardLogger("light_logs", name="mlp_logs_v2")
    early_stop_callback = EarlyStopping(monitor="val_loss", min_delta=1, patience=5, verbose=1, mode="min")
    trainer = pl.Trainer(gpus=0, logger=logger, callbacks=[early_stop_callback])

    dm = MLPDataMod(data_dir=args.data_dir, batch_size=args.batch_size, val_frac=args.val_frac)
    dm.setup()
    data_hparams = dict(
        input_size=dm.input_size,
        batch_size=dm.batch_size,
        val_frac=dm.val_frac
    )

    model = MLP(data_hparams=data_hparams, learn_rate=args.learn_rate, numnodes_hiddenlayers=args.nodes_hlayers,
                dropout_vals=args.drop_vals)

    if args.lr_find:
        lr_finder = trainer.tuner.lr_find(model, dm)
        fig = lr_finder.plot(suggest=True, show=True)
        fig.show()

        model = MLP(data_hparams=data_hparams, learn_rate=lr_finder.suggestion(),
                    numnodes_hiddenlayers=args.nodes_hlayers, dropout_vals=args.drop_vals)

        learn_rate = lr_finder.suggestion()

    trainer.fit(model, dm)
