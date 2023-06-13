import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std

drp = ['structure_id']  # cols to drop from data
s = 0
# Get train data
df = pd.read_csv('./mlp_data_v6.csv').drop(drp, axis=1)
df = df.iloc[:, :df.columns.get_loc('s_bias') + 1]

if s:
    # Std feats
    distinfo = dict()
    for el in df.columns[:-1]:
        mean, std = df[el].mean(), df[el].std()

        if std > 0:
            df[el] = (df[el] - mean) / std
        distinfo[el] = (mean, std)

x, y = df.iloc[:, :df.columns.get_loc('s_bias')], df['s_bias']
sm.add_constant(x)

# Fit data
mod = sm.WLS(y, x, weights=np.array([0.2]*len(df)))
#mod = sm.(y, x, M=sm.robust.norms.HuberT())
res = mod.fit()
print(res.summary())

# Get test data
df = pd.read_csv('./crambin_data.csv').drop(drp, axis=1)
df = df.iloc[:, :df.columns.get_loc('s_bias') + 1]

if s:
    # Std feats
    for el in df.columns[:-1]:
        mean, std = distinfo[el][0], distinfo[el][1]

        if std > 0:
            df[el] = (df[el] - mean) / std

x, y = df.iloc[:, :df.columns.get_loc('s_bias')], df['s_bias']

# Predict
x = sm.add_constant(x)
ypred = res.predict(x)
print(ypred)


# Plotting
import seaborn as sns
sns.set(style='white')

# Plot actual and predicted entropy values
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, dpi=300)
fig.set_size_inches(10, 10)

inds = range(len(y))

sns.scatterplot(x=inds, y=y, ax=ax1, label='Actual')
sns.scatterplot(x=inds, y=ypred, ax=ax2, label='Prediction')
sns.scatterplot(x=inds, y=y, ax=ax3, label='Actual')
sns.scatterplot(x=inds, y=ypred, ax=ax3, label='Prediction')

for ax in (ax1, ax2, ax3):
    ax.legend()

fig.suptitle('Actual & Predicted Configurational Entropy')
fig.supxlabel('Row Number')
fig.supylabel('Entropy')
fig.show()


# Find and plot the residuals for the fit
residuals = [ypred[i] - y[i] for i in inds]
fig, (ax1, ax2) = plt.subplots(2, 1, dpi=300)
fig.set_size_inches(6, 8)

sns.scatterplot(x=inds, y=residuals, ax=ax1)
sns.histplot(data=residuals, ax=ax2)

ax1.set(xlabel='Row number', ylabel='Residual', title='Residual Plot')
ax2.set(xlabel='Residual', ylabel='Count', title='Residual Distribution')

fig.show()

pe = [100 * (ypred[i] - y[i]) / y[i] for i in inds]
fig, (ax1, ax2) = plt.subplots(2, 1, dpi=300)
fig.set_size_inches(6, 8)

sns.scatterplot(x=inds, y=pe, ax=ax1)
sns.histplot(data=pe, ax=ax2)

ax1.set(xlabel='Row number', ylabel='Percent error', title='Percent Error Plot')
ax2.set(xlabel='Percent error', ylabel='Count', title='Percent Error Distribution')

fig.show()


