# -*- coding: utf-8 -*-
"""
filetools
=========

Tools for handling files

by Xu Ruqing, created Jul 2015
"""

import numpy as np
import os.path

def check_create_dir(newdir,basedir=None,makenew=True):
    '''
    Check if newdir exists (relative to basedir);
    create newdir (recursively) if not exist, unless makenew is set to False.
    
    Gives RuntimeWarning if newdir is the same as basedir.
    
    Returns absolute path of newdir
    '''
    import os.path

    ## join newdir to basedir if newdir is not absolute
    if not os.path.isabs(newdir):
        if basedir is None:
            raise RuntimeError("New dir is relative and not base dir provided.")
        elif not os.path.exists(basedir):
            raise RuntimeError("Invalid base directory.")
        else:
            newdir = os.path.join(basedir,newdir)
        
    ## make new dir if it doesn't exist
    if (not os.path.exists(newdir)) and makenew:
        from os import makedirs
        makedirs(newdir)
    elif (basedir is not None) and \
         os.path.abspath(basedir)==os.path.abspath(newdir):
        raise RuntimeWarning('New directory is the same as base directory.')
    
    return newdir
            


class FileSet:
    '''
    A class that defines a set of files, also an iterator

    Usage example:

        inputfiles1 = FileSet('/data','scan_{}.h5',[1,2])
        inputfiles2 = FileSet('/recondata','scan_{}_{}.h5', \
                              (list(range(1,4)),np.array(range(3))))
        for filepath in inputfiles2:
            print(filepath)
        print(inputfiles1[-1])
    '''

    def __init__(self,folder,name_pattern,num_ranges):
        ### input check ###
        assert isinstance(folder,type('a'))
        assert isinstance(name_pattern,type('a'))
        self.folder = str(folder)
        self.namepattern = str(name_pattern)
        __num_of_fields = min(name_pattern.count('{'),name_pattern.count('}'))

        def __check_1d_int_input(p):
            if isinstance(p,np.ndarray):
                # if numpy array, check if it's only 1-d and all data are integers
                assert p.ndim == 1
                assert issubclass(p.dtype.type,np.integer)
                assert set(p>=0) == {True}
            elif isinstance(p,list):
                # if list, check if all elements are integers
                isgood = [isinstance(i,int) and i>=0 for i in p]
                assert len(set(isgood)) == 1
            else:
                raise AssertionError

        # num_ranges is either a 1-d ndarray or a list or a tuple)
        if isinstance(num_ranges,(np.ndarray,list)):
            # if array or list, should be non-negative integers
            __check_1d_int_input(num_ranges)
            self.ranges = (list(num_ranges),)    # store as a tuple of list
        elif isinstance(num_ranges,tuple):
            # if tuple ,should have only 2 elements, each being an array or list
            assert len(num_ranges)==2
            for one_range in num_ranges:
                try:
                    __check_1d_int_input(one_range)
                except:
                    raise AssertionError
            self.ranges = (list(num_ranges[0]),list(num_ranges[1]))

        assert __num_of_fields == len(self.ranges)

        #self.range_sizes = (len(self.ranges[0]),) if self.num_of_fields == 1   \
        #                    else (len(self.ranges[0]),len(self.ranges[1]))
        #self.num_of_files = len(self.ranges[0]) if self.num_of_fields == 1   \
        #                    else len(self.ranges[0])*len(self.ranges[1])

    def get_num_of_fields(self):
        return min(self.namepattern.count('{'),self.namepattern.count('}'))

    def get_range_sizes(self):
        return (len(self.ranges[0]),) if self.num_of_fields == 1   \
                else (len(self.ranges[0]),len(self.ranges[1]))

    def get_num_of_files(self):
        return len(self.ranges[0]) if self.num_of_fields == 1   \
               else len(self.ranges[0])*len(self.ranges[1])

    num_of_fields = property(get_num_of_fields)

    range_sizes = property(get_range_sizes)

    num_of_files = property(get_num_of_files)



    def __iter__(self):
        self.i = 0
        if self.num_of_fields == 2:
            self.j = 0
        return self

    def __next__(self):
        '''return the next file path'''

        ndim2 = self.range_sizes[1] if self.num_of_fields == 2 else 1
        if self.i < self.num_of_files:
            ii = self.i // ndim2
            jj = self.i % ndim2
    
            iii = self.ranges[0][ii]
            jjj = self.ranges[1][jj]  if self.num_of_fields==2 else  0

            self.i += 1
            return os.path.join(self.folder,self.namepattern.format(iii,jjj))
        else:
            raise StopIteration

    def __getitem__(self,i):
        '''returns path of the i-th file'''
        ## check boundary ##
        # allow [-1] style indexing, but cannot exceed total N of files
        if i>=self.num_of_files or i<-self.num_of_files:
            raise IndexError('File index out of range.')

        ## convert index to a number between 0 and N-1 ##
        index = i % self.num_of_files

        ## now find corresponding file ##
        ndim2 = self.range_sizes[1] if self.num_of_fields == 2 else 1
        ii = index // ndim2
        jj = index % ndim2
        iii = self.ranges[0][ii]
        jjj = self.ranges[1][jj]  if self.num_of_fields==2 else  0
        return os.path.join(self.folder,self.namepattern.format(iii,jjj))


if __name__ == '__main__':
    files = FileSet('c','scan_{}.h5',list(range(4)))
    for filename in files:
        print(filename)

    files = FileSet('c','scan_{}_{}.h5',(list(range(2)),np.array(list(range(2)))))
    for filename in files:
        print(filename)
