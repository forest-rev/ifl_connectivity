# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 11:36:12 2019
Create list of urls to give to wget for download
@author: pj276
"""
#%%
import os, glob

#%%
for i in ['af','sa','as']:
    flist = glob.glob('/scratch/pj276/ifl_corridors/gfc/gfc1_0/gfc1_0_' + i + '/' + 'Hansen_GFC2013_datamask_*')
    flist = [j for j in flist if j.endswith('.tif')]
    flist = [os.path.basename(j).split('_')[-2] + '_' + os.path.basename(j).split('_')[-1].replace('.tif','') for j in flist]
    
    # Write list of urls to text file for wget automated download
    with open('/scratch/pj276/ifl_corridors/gfc/downloadlist2_' + i + '.txt', 'w') as afile:
        for k in flist:
            # For some Asia files, I had to change the extent to avoid wrapping.
            # I added a v2 to the file names and these shouldn't be listed here.
            if not 'v2' in k:
                url = 'https://storage.googleapis.com/earthenginepartners-hansen/GFC-2018-v1.6/Hansen_GFC-2018-v1.6_treecover2000_' + k + '.tif'
                afile.write("%s\n" % url)


                


