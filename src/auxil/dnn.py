#!/usr/bin/env python
#******************************************************************************
#  Name:     dnn.py
#  Purpose:  object class for supervised image classification with  NN (tensorflow)
#  Usage:    
#     from dnn import Dnn
#
# (c) Mort Canty 2019

import numpy as np  
import pandas as pd
import matplotlib.pyplot as plt
import tensorflow as tf 
from tensorflow.keras import layers

#tf.logging.set_verbosity('ERROR')   
    
class Dnn(object):    
    '''High-level TensorFlow (keras) Dnn classifier'''
    def __init__(self,Ls=[10],n_classes=10,learning_rate=0.01):
#      setup the network architecture  
        self._dnn = tf.keras.models.Sequential()
        self._history = None
#      hidden layers        
        for L in Ls:
            self._dnn \
             .add(layers.Dense(L,activation='relu'))
#      output layer
        self._dnn \
          .add(layers.Dense(n_classes,
                            activation='softmax'))       
#      initialize                             
        self._dnn.compile(
                optimizer=tf.keras.optimizers.SGD(learning_rate=learning_rate),
                metrics=['accuracy'],
                loss='categorical_crossentropy')
        
    def train(self,Gs,ls,epochs=10):
        n_split = (2*ls.shape[0])//3
        self._Gs_train = Gs[:n_split,:]
        self._Gs_valid = Gs[n_split:,:]
        self._ls_train = ls[:n_split,:]
        self._ls_valid = ls[n_split:,:]
        try:           
            self._history = self._dnn.fit(self._Gs_train,self._ls_train,
                        epochs=epochs,verbose=0,
                        validation_data=(self._Gs_valid,self._ls_valid))
            return True 
        except Exception as e:
            print( 'Error: %s'%e ) 
            return None        
        
    def history(self,sfn=None):
        pd.DataFrame(self._history.history).plot(figsize=(8,5))
        plt.grid(True)
        plt.gca().set_ylim(0,1)
        if sfn is not None:
            plt.savefig(sfn,bbox_inches='tight')    
        plt.show()                                               
        
    def classify(self,Gs):     
#      predict new data                       
        Ms = self._dnn.predict(Gs)
        cls = np.argmax(Ms,1)+1 
        return (cls,Ms)

    def test(self,Gs,ls):
        m = np.shape(Gs)[0]
        classes, _ = self.classify(Gs)
        classes = np.asarray(classes,np.int16)
        labels = np.argmax(np.transpose(ls),axis=0)+1
        misscls = np.where(classes-labels)[0]
        return len(misscls)/float(m)  
    
    def save(self,path):
        self._dnn.save(path)
        
    def restore(self,path):
        self._dnn = tf.keras.models.load_model(path)
    
if __name__ == '__main__':
#  test on random data    
    Gs = 2*np.random.random((1000,3)) -1.0
    ls = np.zeros((1000,4))
    for l in ls:
        l[np.random.randint(0,4)]=1.0 
    cl = Dnn(Gs,ls,[10,10])  
    if cl.train():
        classes, probabilities = cl.classify(Gs)
        print( classes )