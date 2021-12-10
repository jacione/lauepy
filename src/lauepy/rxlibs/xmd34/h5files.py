# -*- coding: utf-8 -*-
"""
h5files
=======
Tools for handling HDF5 image files taken at 34ID-E

Includes:
---------
* class XmdH5Data
* class WireStepScanFileSet

By Ruqing Xu, created Aug.2015
"""

import numpy as np
import h5py
import os.path

class XmdH5Data(object):
    '''
    Class that loads data from an HDF5 image file collected at beamline 34ID-E.

    Usage
    -----
    ::
  
        object = XmdH5Data(h5filepath, **kwargs)
        
    Parameters
    ^^^^^^^^^^
        h5filepath : string
            string that contains full path of the data file, must have ".h5" extention.

    Accepted keyword arguments (kwargs):
    
        load_metadata : boolean
            Whether to load metadata such as sample positions, beam mode etc. 
            If True, metadata will be stored as attributes, see Attributes/Properties 
            section below for details. Default = True
        load_main_array : boolean
            Whether to load and keep the full data array in memory;
            this option is useful when dealing with large 3D dataset from fly scans.
            Default is False.
        readonly : boolean
            Whether to open the file read-only;
            If true, some writing methods will return Exceptions if called.
            Default is True.
        keep_file_open: boolean
            Whether to keep the file object open after data is loaded;
            If true, self.filehandle will be an open h5py.File() instance that
            the user can use to custom-handle the data file. In this case the user
            is strongly suggested to use the "with" construct (see Usage Suggestions
            below).
            Default is False
    Usage suggestions
    ^^^^^^^^^^^^^^^^^
        If keep_file_open is True, then the file object needs to be properly closed in the end.
        The best practice is to use this class in a "with" construct, which will automatically
        do the file closing upon finishing its contents, for example::

                with XmdH5Data('/data/test/a.h5',keep_file_open=True) as image_a:
                    print numpy.amax(image_a.data)
                    ...

        If not used in a "with" construct, remember to use the .closeout() method, which will
        close the file and remove the .data array from memory.

    Other Notes
    ^^^^^^^^^^^
        * Main data array is loaded into memory as 8-byte float numbers instead of integers.
          The methods that writes data back to the file, however, writes in the same
          integer format as in the original file.
        * Uses the **h5py** package to open .h5 files

    Attributes / Properties
    -----------------------
        (all arrays are numpy arrays)
        
        .filepath
            original path of the file
        .filename
            (property,read-only) file name part of .filepath
        .filedir      
            (property,read-only) directory part of .filepath
        .filehandle
            h5py.File() object if file is opened
        .file_is_open
            (property,read-only) boolean, whether h5 file is open via self.filehandle
        .load_main_array:
            (property) boolean, whether full data array is stored in memory or not.
            if value is changed by user, the data array will be automatically loaded
            or removed upon the change.
        .data: 
            (property) 2D or 3D array of 8-byte float type containing the main data;
            =None if store_full_data=False
        .scanmode
            (property,read_only) 'fly' or 'step'
        .data_dims
            (property,read-only) tuple that gives the shape of the main data array
        .dimX, .dimY
            (property,read-only) image dimensions
        .num_of_frames
            (property,read-only) number of frames for fly scans, = 1 for step scans.
            
        (The following won't exist if load_metadata=False upon initialization)
        
        .I0
            (property) alias of self.I_middle
        .I_middle, .I_final, .I_start
            beamline scaler counts
        .beammode
            (property,read-only) string, 'monochromatic' or 'white' etc.
        .energy
            beam energy in keV
        .depth
            depth of reconstruction, only if present in the hdf5 file
        .sampleX, .sampleY, .sampleZ
            sample positions
        .sampleH, .sampleF
            (property,read_only) number computed from sampleY and sampleZ
        .wirebaseX, .wirebaseY, .wirebaseZ
            wire base motor positions, numbers only, unlike .wireX/Y/Z
        .wireX, .wireY, .wireZ
            wire positions, numbers for step scans, 1D arrays for fly scans
        .wirescanPos
            position of wire scan motor (usually H), number or 1D array
        .wireH, .wireF
            (property,read_only) numbers or 1D arrays computed from wireY and wireZ
        .roi_startX, .roi_startY
            ROI starting indices
        .roi_endX, .roi_endY
            (property,read_only) computed from roi_start indices and image dimensions
        .detectorID
            (property,read-only) string, ID of detector
        .exposure
            exposure time in seconds

    Methods
    --------
        .loaddata()
            opens file and reads data fields
        .reload()
            reload all fields from file
        .closeout()
            close file object, delete .data array and set store_full_data to False;
            object still remains with other attributes intact.
        .get_full_data()
            returns the 2D or 3D data array
        .get_frame(n)
            returns the n-th frame from fly scan data set, n starts from 0
        .write_data(newdata=None)
            writes data array to file, if 'newdata' is provided, write it;
            if not, write self.data .
            Raises error if file mode is read-only, or if data is not in memory.
        .write_frame(n,newframe)
            writes 2D array newframe to the n-th frame in the file.
            Raises error if file mode is read-only
        .copyfile(newfilepath)
            copy current image file to a new path, newfilepath can be just
            another folder or can also include a new file name
        .process_data(proc_func,*fargs,**fkwargs)
            generic function to allow user to call a custom fucntion to 
            process the data array, proc_func should take a 2D numpy array
            as 1st argument
        .process(proc_func,*fargs,**fkwargs)
            generic function to allow user to call a custom fucntion to 
            process image, proc_func should take a XmdH5Data instance
            as the 1st argument; this is useful if the process function
            also needs image metadata besides the main data array
    '''
    def __init__(self,h5path,load_metadata=True,load_main_array=None, \
                      readonly=True,keep_file_open=False,**other_args):
        assert isinstance(h5path,type('a'))
        assert h5path[-3:] == '.h5'  # extention must be '.h5'

        self.filepath = h5path
        self._load_metadata = load_metadata
        if load_main_array is not None:
            self._store_full_data = load_main_array
        elif 'store_full_data' in other_args: # for backwd compatibility
            self._store_full_data = other_args['store_full_data']
        else:
            self._store_full_data = False  # default value
        
        self._keep_file_open = keep_file_open

        self._filemode = 'r' if readonly else 'r+'
        self.loaddata()

    def loaddata(self):
        '''
        Opens h5 file and reads data fields.

        File is closed in the end if _keep_file_open is False.
        '''

        try:
            self.filehandle = h5py.File(self.filepath,mode=self._filemode)
        except:
            raise OSError('Problem openning file {}'.format(self.filepath))

        ### read data & mode ###
        self._data_type = self.filehandle['entry1']['data']['data'].dtype

        # Note: data matrix stored in hdf5 has its [i,j] element corresponding 
        # to the i-th row and j-th column. While the X-direction is defined to
        # be in the "horizontal" direction of the matrix and thus corresponds 
        # to increase of the column index, j, and hence Y-direction to row  
        # index i. 
        # to allow intuitive handling of data elements, i.e., let X be the 1st
        # axis and Y be the 2nd axis and refer to data elements by data[ix,iy],
        # we impose a transpose of the last 2 dimensions of the data array from
        # the hdf5 file.
        data_shape = list(self.filehandle['entry1']['data']['data'].shape)
        data_shape[-2],data_shape[-1] = data_shape[-1],data_shape[-2]
        self._data_dims = data_shape
        dims = list(range(len(self._data_dims)))
        dims[-2],dims[-1] = dims[-1],dims[-2]
        if self._store_full_data:
            self._data = np.array(self.filehandle['entry1']['data']['data'],  \
                                  dtype='f8').transpose(dims)

        ### read metadata ###
        if self._load_metadata:
            try:
                self._beammode = self.filehandle['entry1']['microDiffraction'
                                                 ]['MonoMode'][0]
            except:
                # older h5 files does not have this entry
                self._beammode = None
    
            ## sample positions and beam energy ##
            self.energy = self.filehandle['entry1']['sample']['incident_energy'][0]
            self.sampleX = self.filehandle['entry1']['sample']['sampleX'][0]
            self.sampleY = self.filehandle['entry1']['sample']['sampleY'][0]
            self.sampleZ = self.filehandle['entry1']['sample']['sampleZ'][0]
    
            ## wire positions ##
            if self.filehandle['entry1'].__contains__('wire'):
                _wire = self.filehandle['entry1']['wire']
                self.wirebaseX = _wire['wirebaseX'][0]
                self.wirebaseY = _wire['wirebaseY'][0]
                self.wirebaseZ = _wire['wirebaseZ'][0]
                if len(self._data_dims)==2:  # is step scan
                    self.wireX = _wire['wireX'][0]
                    self.wireY = _wire['wireY'][0]
                    self.wireZ = _wire['wireZ'][0]
                    self.wirescanPos = _wire['wirescan'][0]
                elif len(self._data_dims)==3: # is fly scan
                    self.wireX = np.array(_wire['wireX'],dtype=_wire['wireX'].dtype)
                    self.wireY = np.array(_wire['wireY'],dtype=_wire['wireY'].dtype)
                    self.wireZ = np.array(_wire['wireZ'],dtype=_wire['wireZ'].dtype)
                    self.wirescanPos = np.array(_wire['wirescan'],           \
                                              dtype=_wire['wirescan'].dtype)
                    assert len(self.wireX) == self._data_dims[0]+2
    
            ## other info ##
            self._detectorID = self.filehandle['entry1']['detector']['ID'][0]
            self.roi_startX = self.filehandle['entry1']['detector']['startx'][0]
            self.roi_startY = self.filehandle['entry1']['detector']['starty'][0]
            self.exposure = self.filehandle['entry1']['detector']['exposure'][0]
            self.I_middle = self.filehandle['entry1']['monitor']['I0_calc'][0]
            self.I_final = self.filehandle['entry1']['monitor']['I_final_calc'][0]
            self.I_start = self.filehandle['entry1']['monitor']['I_start_calc'][0]
            if self.filehandle['entry1'].__contains__('depth'):
                self.depth = self.filehandle['entry1']['depth'][0]

            ### end of reading metadata
                
        if not self._keep_file_open:
            self.filehandle.close()

    def reload(self):
        '''
        Reload data from file.
        If file is open, close and re-open it.
        '''
        # close & reopen file
        if self.file_is_open:
            self.filehandle.close()
        self.loaddata()

    @property
    def filename(self):
        return os.path.basename(self.filepath)
        
    @property
    def filedir(self):
        return os.path.dirname(self.filepath)

    @property
    def file_is_open(self):
        return hasattr(self,'filehandle') and self.filehandle.__contains__('entry1')

    @property
    def load_main_array(self):
        '''
        True or False: whether the full data array is stored in memory;

        if set to False, the "data" attribute will be empty (None);
        if changed to True from False, data array will be automatically loaded from file.
        '''
        return self._store_full_data

    @load_main_array.setter
    def store_full_data(self,value):
        if self._store_full_data == value:
            pass
        elif value==False:
            self._store_full_data = value
            self._data = None
        elif value==True:
            self._store_full_data = value
            if not self.file_is_open:
                self.filehandle = h5py.File(self.filepath,mode=self._filemode)
            dims = list(range(len(self._data_dims)))
            dims[-2],dims[-1] = dims[-1],dims[-2]
            self._data = np.array(self.filehandle['entry1']['data']['data'], \
                                  dtype='f8').transpose(dims)
            if not self._keep_file_open:
                self.filehandle.close()

    @property
    def data(self):
        if self._store_full_data:
            return self._data
        else:
            return None

    @data.setter
    def data(self,newdata):
        if self._filemode == 'r':
            raise IOError('Data set was initialized as read-only!')
        if self._store_full_data:
            if newdata.shape==self.data_dims:
                self._data = newdata
            else:
                raise ValueError('New data array has wrong shape')
        else:
            pass

    @property
    def I0(self):
        if self._load_metadata:
            return self.I_middle
        else:
            return None

    @property
    def beammode(self):
        if self._beammode:
            return self._beammode.decode('utf-8')  #<<< for Py3 only <<<
        else:
            return 'unknown'

    @property
    def scanmode(self):
        if len(self.data_dims)==3:
            return 'fly'
        elif len(self.data_dims)==2:
            return 'step'
        else:
            return 'unknown'

    @property
    def sampleH(self):
        if self._load_metadata:
            return (self.sampleY+self.sampleZ)/np.sqrt(2.0)
        else:
            return None

    @property
    def sampleF(self):
        if self._load_metadata:
            return (self.sampleZ-self.sampleY)/np.sqrt(2.0)
        else:
            return None

    @property
    def wireH(self):
        if self._load_metadata and hasattr(self,'wireZ'):
            return (self.wireY+self.wireZ)/np.sqrt(2.0)
        else:
            return None

    @property
    def wireF(self):
        if self._load_metadata and hasattr(self,'wireZ'):
            return (self.wireZ-self.wireY)/np.sqrt(2.0)
        else:
            return None

    @property
    def roi_endX(self):
        if self._load_metadata:
            return self.roi_startX + self.data_dims[-2] - 1
        else:
            return None

    @property
    def roi_endY(self):
        if self._load_metadata:
            return self.roi_startY + self.data_dims[-1] - 1
        else:
            return None

    @property
    def data_dims(self):
        return self._data_dims

    @property
    def num_of_frames(self):
        if len(self.data_dims)==3:
            return self.data_dims[0]
        elif len(self.data_dims)==2:
            return 1
        else:
            return None

    @property
    def dimX(self):
        return self.data_dims[-2]

    @property
    def dimY(self):
        return self.data_dims[-1]

    @property
    def detectorID(self):
        if self._load_metadata:
            return self._detectorID.decode('utf-8')  #<<< Py3 only <<<
        else:
            return None

    def get_full_data(self):
        '''return full 2d or 3d data array'''
        if self._store_full_data:
            return self._data
        else:
            if not self.file_is_open:
                self.filehandle = h5py.File(self.filepath,mode=self._filemode)
            dims = list(range(len(self._data_dims)))
            dims[-2],dims[-1] = dims[-1],dims[-2]
            _data = np.array(self.filehandle['entry1']['data']['data'],
                             dtype='f8').transpose(dims)
            if not self._keep_file_open:
                self.filehandle.close()
            return _data

    def get_frame(self,frame_No):
        ''' return 1 image frame from a 3d array; frame_No is 0-based '''
        assert isinstance(frame_No,int)
        if frame_No >= self.num_of_frames:
            return None
        elif self.scanmode=='step':
            return self.get_full_data()
        elif self.scanmode=='fly':
            if not self.file_is_open:
                self.filehandle = h5py.File(self.filepath,mode=self._filemode)
            _data = np.array(self.filehandle['entry1']['data']['data'][frame_No],
                             dtype='f8').transpose()
            if not self._keep_file_open:
                self.filehandle.close()
            return _data

    def write_data(self,newdata=None):
        '''
        write data to file 

        if newdata is not provided, self.store_full_data should be True and
        new data should already be stored as self.data
        '''
        if self._filemode == 'r':
            raise IOError('Data set was initialized as read-only!')

        ## if newdata is not provided
        if newdata is None:
            ## use stored data, if not stored, raise error
            _data = self.data
            if _data is None:
                raise ValueError('No data to write!')
        ## if new data is provided
        else:
            ## use new data, check dimensions first
            assert type(newdata) == np.ndarray
            newshape = list(newdata.shape)
            newshape[-2],newshape[-1] = newshape[-1],newshape[-2]
            if newshape != self.data_dims:
                raise ValueError('Data to be written to hdf file has wrong shape')
            _data = newdata

        ## make sure limits are not exceeded
        if issubclass(self._data_type.type, np.integer):
            _limits = np.iinfo(self._data_type)
            _data[_data>_limits.max] = _limits.max
            _data[_data<_limits.min] = _limits.min

        ## open file & write to it, close after done if desired
        if not self.file_is_open:
            self.filehandle = h5py.File(self.filepath,mode='r+')
        # transpose data matrix to match the correct dimensions in hdf5
        dims = list(range(len(self._data_dims)))
        dims[-2],dims[-1] = dims[-1],dims[-2]
        self.filehandle['entry1']['data']['data'][...] = \
                       _data.astype(self._data_type).transpose(dims)
        if self._keep_file_open:
            self.filehandle.flush()
        else:
            self.filehandle.close()


    def write_frame(self,frame_No,newframe):
        '''write a 2D frame to a fly scan file'''
        if self._filemode == 'r':
            raise IOError('Data set was initialized as read-only!')
        if self.scanmode != 'fly':
            raise Exception('Method write_frame() can only be used for flyscan datasets.')

        ## check data shape ##
        assert isinstance(frame_No,int)
        assert type(newframe) == np.ndarray
        if newframe.shape[::-1] != self.data_dims[-2:]:
            raise ValueError('Frame data to be written to hdf file has wrong shape')

        ## check limits ##
        if issubclass(self._data_type.type, np.integer):
            _limits = np.iinfo(self._data_type)
            newframe[newframe>_limits.max] = _limits.max
            newframe[newframe<_limits.min] = _limits.min

        ## open file & write to it, close after done if desired
        if not self.file_is_open:
            self.filehandle = h5py.File(self.filepath,mode='r+')
        _data = self.filehandle['entry1']['data']['data']
        _data[frame_No,:,:] = newframe.astype(self._data_type).transpose()
        if self._keep_file_open:
            self.filehandle.flush()
        else:
            self.filehandle.close()

    def __enter__(self):
        return self

    def __exit__(self,errortype,errorvalue,traceback):
        self.closeout()

    def closeout(self):
        '''
        Close the file and the delete the "self.data" array. The object (self) 
        will continue to exist with other attributes which don't cost much memory.

        This method should be called if object is not used in a "with" 
        construction.
        '''
        if self._keep_file_open:
            self.filehandle.close()
        self.store_full_data = False  # this will also set self._data to None

    def copyfile(self,newfilepath):
        ''' 
        Copy image file to a new file (will overwrite if new file exists)
        '''
        from shutil import copy2
        copy2(self.filepath,newfilepath)
        
    def process(self,proc_func,*fargs,**fkwargs):
        ''' 
        generic function to call user-provided function to do image processing;
        proc_func is user-supplied function that accept an XmdH5Data object as 
        1st argument; *fargs and **fkwargs will be passed to proc_func.
        '''
        assert callable(proc_func)
        proc_func(self,*fargs,**fkwargs)
        
    def process_data(self,proc_func,*fargs,**fkwargs):
        ''' 
        generic function to call user-provided function to do processing over
        the data array only; 
        proc_func is user-supplied function that accept a numpy array as 
        1st argument; *fargs and **fkwargs will be passed to proc_func.
        '''
        assert callable(proc_func)
        if self._data is None:
            raise ValueError('No data to process!')
        else:
            self._data = proc_func(self._data,*fargs,**fkwargs)
            



#from ..misc.filetools import FileSet as __FileSet
from .extutils import filetools as __filetools
__FileSet = __filetools.FileSet

class WireStepScanFileSet(__FileSet):
    '''
    Defines a set of files from a wire step scan.
    
    A subclass of filetools.FileSet, same initialization but with these additional
    attributes/methods:
    
    .check_wire_travels():
        Go through entire file set and find out the file numbers that correspond
        to the beginning of each wire travel.
    .wirestart_indices:
        indices of the files for each beginning of wire travel; crated after 
        calling self.check_wire_travels(). **Note**: this is not the same as 
        the numbers in the filename string, but rather an index counted from 
        0 to number_of_files within this file set.
    .num_wiretravels:
        (property, read-only) the total number of wire travels within this file
        set; if self.wirestart_indices does not exist, then
        self.check_wire_travels() is called automatically.
    '''
    def check_wire_travels(self):
        nfiles = self.num_of_files
        assert nfiles>2
        filepath = self.__getitem__(0)
        with XmdH5Data(filepath) as image1:
            assert image1.scanmode == 'step'
            wirePos0 = image1.wirescanPos
        with XmdH5Data(self.__getitem__(1)) as image1:
            wirePos1 = image1.wirescanPos
            if abs(wirePos1-wirePos0) < 0.25:
                raise Exception('Wire movement between first 2 images < 0.25um.')
                
        wirePos_last = wirePos1
        self.wirestart_indices = [0]
        for i in range(2,nfiles):
            with XmdH5Data(self.__getitem__(i)) as image1:
                wirePos = image1.wirescanPos
            wiredist_last = max(wirePos_last - wirePos0, 0.1)
            if((wiredist_last > 0.25) and \
               ((wirePos - wirePos_last)/wiredist_last < -0.5 )):
                # if wire has moved at all, and has gone backwards with a big step,
                # then record this step as beginning of a new wire scan
                self.wirestart_indices.append(i)
            wirePos_last = wirePos
        
    @property
    def num_wiretravels(self):
        if not hasattr(self,'wirestart_indices'):
            self.check_wire_travels()
        return len(self.wirestart_indices)



        
if __name__ == '__main__':
    filepath = r'D:\temp\2 min\recon\EWscanfly_sort_1_35.h5'
    
    img = XmdH5Data(filepath,store_full_data=True,readonly=False)
    