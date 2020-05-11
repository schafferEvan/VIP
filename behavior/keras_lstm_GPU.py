import numpy as np
import pandas as pd 
import tensorflow as tf
import sys
from keras.layers import Dense, Dropout
from keras.models import Sequential,model_from_yaml
import keras 
from numpy import mean, std, dstack
from keras.layers import Flatten, LSTM
from keras.utils import to_categorical
from sklearn.metrics import classification_report

#/mnt/ceph/users/evanschaffer/data/classifier/
# load data from the csv files
print('\nLoading files..\n')
train_csv = pd.read_csv('/mnt/ceph/users/evanschaffer/data/classifier2p0/train3.csv')
# train_csv = pd.read_csv('/Volumes/ESSbehavior/classifier2p0/train3.csv')
train_csv = train_csv.values

ground_csv = pd.read_csv('/mnt/ceph/users/evanschaffer/data/classifier2p0/groundtruth3.csv')
# ground_csv = pd.read_csv('/Volumes/ESSbehavior/classifier2p0/groundtruth3.csv')
ground_csv = ground_csv.values

# separate into input x and target y dataframes...
x_all = train_csv[:,0:485]
y_all = train_csv[:,485]
gx_all = ground_csv[:,0:485]
gy_all = ground_csv[:,485]

# format testing ratios
TEST_SIZE = 0
x_test = x_all[:TEST_SIZE,:]
y_test = y_all[:TEST_SIZE]
x_train = x_all[TEST_SIZE:,:]
y_train = y_all[TEST_SIZE:]


#normalization of data
m1 = x_all.mean()
for i in x_all:
    i -= m1
m2 = gx_all.mean()
for i in gx_all:
    i -= m1

x_all = x_all/x_all.max()
gx_all = gx_all/gx_all.max()

# # scale data
# USE_SCALER = True
# if USE_SCALER == True:
#     from sklearn.preprocessing import StandardScaler
#     scaler = StandardScaler()
#     scaler.fit(x_train)   # Fit only to the training dataframe
#     x_train = scaler.transform(x_train)

# convert data to 3D array
print(x_train.shape[1])
x_train = x_train.reshape((x_train.shape[0],x_train.shape[1],1))
x_test = x_test.reshape((x_test.shape[0],x_test.shape[1],1))
y_train = to_categorical(y_train)
gx_all = gx_all.reshape((gx_all.shape[0],gx_all.shape[1],1))
gy_all = to_categorical(gy_all)

model = Sequential()
model.add(LSTM(100, input_shape=(x_train.shape[1],x_train.shape[2])))
model.add(Dropout(0.5))
model.add(Dense(100, activation='relu'))
model.add(Dense(2, activation='sigmoid'))
model.compile(loss='binary_crossentropy', optimizer='rmsprop', metrics=['accuracy'])

# protect against early overfitting by reducing the learning rate.
reduce_lr = keras.callbacks.ReduceLROnPlateau(monitor='loss', factor=0.2,patience=3, min_lr=0.001)

# fit the model to the data
print('\nFitting the model..\n')
model.fit(x_train, y_train, epochs=15, batch_size=64,callbacks=[reduce_lr])

# serialize model to YAML
model_yaml = model.to_yaml()
with open("model.yaml", "w") as yaml_file:
    yaml_file.write(model_yaml)
# serialize weights to HDF5
model.save_weights("model.h5")
print("Saved model to disk")

# score the model
print('\nScoring the model..\n')
# score = model.evaluate(x_test, y_test, batch_size=128,)
x_all = x_all.reshape((x_all.shape[0],x_all.shape[1],1))
y_all = to_categorical(y_all)
score = model.evaluate(x_all, y_all, batch_size=128,)
print(model.metrics_names)
print(score)

#print(classification_report(gy_all,model.predict_classes(gx_all)))


# Stage 2 ------------------------------------------------------------------------------------------------------------

# load YAML and create model
yaml_file = open('model.yaml', 'r')
loaded_model_yaml = yaml_file.read()
yaml_file.close()
model = model_from_yaml(loaded_model_yaml)

# load weights into new model
model.load_weights("model.h5")
print("Loaded model from disk")


# unknown_csv = unknown_csv.reshape((unknown_csv.shape[0],unknown_csv.shape[1],1))
x = model.predict_classes(gx_all)
np.savetxt('all_predictions_e15gpu.csv', x, fmt='%d')
