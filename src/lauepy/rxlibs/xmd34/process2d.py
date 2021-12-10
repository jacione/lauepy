# -*- coding: utf-8 -*-
"""
image processing tools
"""

import numpy as np

class ImageEnhancer:
    '''
    A class that defines an image-processing sequence, which includes:
        1. do a median filter (with custom size) to get a rough background
        2. do a difference-check between the background and signal, on each pixel
           calculate the ratio of (itensity - background)/sqrt(background);
           generate a mask based on a cutoff value of this ratio
           mask = 1 for pixels that belog to a peak
        3. expand the mask size by a given number of pixels to include more tails
        4. (optional) refine the background by doing interplation at the mask=1
           pixels from the background value of the mask=0 pixels
    Usage:
      * Initialize the processor:
        processor = ImageEnhancer(median_filter_size=5,diff_check_cutoffratio=4.0,
                                  pixels_expand=1,ignore_singular=True, subtract_bgd=True)
      * run the processor:
        result = processor.run(dataarray)
    '''
    def __init__(self,median_filter_size=5,diff_check_cutoffratio=4.0, \
                 pixels_expand = 1, ignore_singular=True, subtract_bgd=False):
        self.median_filter_size = median_filter_size
        self.diff_check_cutoffratio = diff_check_cutoffratio
        self.ignore_singular = ignore_singular
        self.subtract_bgd = subtract_bgd
        
        self.kernel_single_next_px = np.array([[0,1,0],[1,0,1],[0,1,0]],dtype='f4')
        self.kernel_expand = self.__gen_expand_kernel(pixels_expand)
        
    def gen_mask_n_bgd(self,data,calc_bgd=True):
        '''
        generate mask & background by using median-filter followed by 
        difference-check, and finally a 2d interplation
        '''
        from scipy.ndimage.filters import median_filter as __median_filter
        from scipy.ndimage import convolve as __convolve
        assert isinstance(data,np.ndarray)
        bgd1 = __median_filter(data,self.median_filter_size)

        ## look for pixels that are too different from median_filtered values ##
        diff_mask = data - bgd1 > np.sqrt(bgd1)*self.diff_check_cutoffratio

        ## ignore if a single pixel stands out from background ##
        # 1. convolve mask with a window function
        tempmask = __convolve(diff_mask,self.kernel_single_next_px,mode='mirror') > 1.0e-5
        # 2. do an AND with original mask
        tempmask = np.logical_and(tempmask,diff_mask)

        ## expand the mask by given number of pixels ##
        self.mask = __convolve(tempmask,self.kernel_expand,mode='mirror') > 1.0e-5
        
        ## if needed, calc background at masked pixels with interpolation
        if calc_bgd:
            return self._intp_masked(bgd1)

    def _mask_is_empty(self):
        '''test if mask is empty'''
        assert hasattr(self,'mask')
        return self.mask.sum()==0
    
    def _intp_masked(self,data):
        '''
        do interplation to get pixels within a mask based on pixels outside of mask
        '''
        assert hasattr(self,'mask')
        from scipy.interpolate import griddata as __grid_intp

        if self._mask_is_empty():
            ## do nothing if mask is empty ##
            return data
        else:
            ## do interpolation in general ##
            mask_n = np.logical_not(self.mask)
            ydim,xdim = data.shape
            xx,yy = np.meshgrid(np.arange(ydim,dtype='f4'),np.arange(xdim,dtype='f4'), \
                               indexing='ij')  #<<< be extra ware of the "indexing" KW <<<
            new = __grid_intp((xx[mask_n],yy[mask_n]),data[mask_n],(xx,yy),method='linear')
            return new
        
    def gen_bgd_withmask(self,data):
        '''
        generate background by using median-filter followed by 
        a 2d interplation, with an existing mask array
        '''
        assert hasattr(self,'mask')
        assert isinstance(data,np.ndarray)
        assert self.mask.shape == data.shape
        from scipy.ndimage.filters import median_filter as __median_filter
        bgd1 = __median_filter(data,self.median_filter_size)
        
        ## calc background with interpolation
        return self._intp_masked(bgd1)
    
    def run(self,data,use_existing_mask=False,output_bgd=False):
        '''
        Apply the full procedure to a data array
        '''
        if use_existing_mask:
            ## sanity check ##
            assert self.mask.shape == data.shape
            
        #### subtract background, or just return the masked pixels ####
        if self.subtract_bgd:
            if not use_existing_mask:
                bgd = self.gen_mask_n_bgd(data)
            else:
                bgd = self.gen_bgd_withmask(data)
            
            result = np.ones(data.shape,dtype=data.dtype)
            result[self.mask] = data[self.mask] - bgd[self.mask]
            result[result<1] = 1  # all-zero arrays may cause problem in LaueGo
            
            result = (result,bgd) if output_bgd else result
                
        else:
            ### if bgd is not needed to be subtracted from data ###
            if not use_existing_mask:
                bgd = self.gen_mask_n_bgd(data,calc_bgd=output_bgd)
            else:
                bgd = self.gen_bgd_withmask(data)
                
            if self._mask_is_empty():
                # do nothing if mask is empty
                result = (data,bgd) if output_bgd else data
            else:
                ## return masked array, +background if needed ##
                result = np.ones(data.shape,dtype=data.dtype)
                result[self.mask] = data[self.mask]
                result = (result,bgd) if output_bgd else result

        return result

    @staticmethod
    def __gen_expand_kernel(R):
        '''
        generate a kernel matrix that is used to convolve with an existing 2d data matrix to 
        expand a certain mask by radius R
        '''
        sizeh = int(np.floor(R))
        size =  sizeh * 2 + 1
        kernel = np.zeros((size,size))
        for i in range(size):
            for j in range(size):
                r = np.sqrt((i - sizeh)**2 + (j-sizeh)**2)
                kernel[i,j] = 1 if r < R else 0
        return kernel