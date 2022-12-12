# run
# python cleanup.py

# this pyhton script delete the .dSYM debug folder and all .exe/.ipynb_checkpoints files under the current directory

import os
import shutil


def cleanup(parent_dir = None):
    if (parent_dir == None):
        curr_ls = os.listdir()
    else:
        curr_ls = os.listdir(parent_dir)
    
    for fileName in curr_ls:
        if (parent_dir == None):
            relative_dir = fileName
        else:
            relative_dir = parent_dir + '/' + fileName
        
        if (os.path.isfile(relative_dir) & 
            ((".exe" in fileName) | 
             ("tempCodeRunnerFile" in fileName))):
            
            os.remove(relative_dir)
            continue
        
        if (os.path.isdir(relative_dir) & 
            ((".dSYM" in fileName) | 
             (".ipynb_checkpoints" in fileName))):
            shutil.rmtree(relative_dir, ignore_errors=True)
            continue
        
        if (os.path.isdir(relative_dir)):
            cleanup(relative_dir)
            
    return


cleanup()