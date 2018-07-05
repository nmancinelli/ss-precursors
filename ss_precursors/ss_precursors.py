#!/bin/env python
import pickle

class Stack:
    def read(filename):
        fin = open(filename,'rb')
        stk = pickle.load(fin)
        return stk
    
    def __init__(self):
        self.amplitude = []
        self.standard_error = []
        self.depths = []
        self.metadata = {}

    def save(self,filename):
        fout = open(filename,'wb')
        pickle.dump(self, fout)

    def print(self):
        print('Amps: ',self.amplitude)
        print('SE  : ',self.standard_error)
        print('MD  : ',self.metadata)



def test():
    t = Stack()
    t.amplitude = [1, 2, 3]
    t.standard_error = [4, 5, 6]
    t.metadata = {1: 2, 3: 4}
    t.save('test.pickle')

    v  = Stack().read('test.pickle').print()


if __name__ == "__main__":
    test()