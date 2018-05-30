#!/bin/env python
import pickle

class Stack:
    def __init__(self):
        self.amplitude = []
        self.standard_error = []
        self.metadata = {}

    def save(self,filename):
        fout = open(filename,'wb')
        pickle.dump(self, fout)

    def read(self, filename):
        fin = open(filename,'rb')
        return pickle.load(fin)

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