#!/usr/bin/env python
#******************************************************************************
#  Name:     dnn.py
#  Purpose:  object class for supervised image classification with  NN (tensorflow)
#  Usage:    
#     from dnn import Dnn
#
# (c) Mort Canty 2019

import numpy as np  
import tensorflow as tf 
from tensorflow.keras import layers

#tf.logging.set_verbosity('ERROR')   
    
class Dnn(object):    
    '''High-level TensorFlow (keras) Dnn classifier'''
    def __init__(self,Gs,ls,Ls,epochs=100,learning_rate=0.01):
#      setup the network architecture  
        self._Gs = Gs   
        n_classes = ls.shape[1]
        self._labels = ls
        self._epochs = epochs
        self._dnn = tf.keras.Sequential()
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
                optimizer=tf.optimizers.SGD(learning_rate=learning_rate),
                loss='categorical_crossentropy')
        
    def train(self):
        try:           
            self._dnn.fit(self._Gs,self._labels,
                        epochs=self._epochs,verbose=0)
            return True 
        except Exception as e:
            print( 'Error: %s'%e ) 
            return None             
        
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