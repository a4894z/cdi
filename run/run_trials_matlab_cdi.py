import os
import subprocess
import sys
import shutil

# os.rename("path/to/current/file.foo", "path/to/new/destination/for/file.foo")
# shutil.move("path/to/current/file.foo", "path/to/new/destination/for/file.foo")
# os.replace("path/to/current/file.foo", "path/to/new/destination/for/file.foo"

#import os
#import shutil as sh
#
#root_path = "./folder_A"
#dest_path = "./folder_B"
#
#for dirpath, dnames, fnames in os.walk(root_path):    
#    for f in fnames:
#        if f.endswith(".jpg"):
#            source_file_path =  os.path.join(dirpath, f)
#            dest_file_path   =  os.path.join(dest_path, f)
#            sh.copyfile(source_file_path, dest_file_path)

def main( ):

    #==================================================================================================
    
    os.system( 'export LD_LIBRARY_PATH=~/Documents/MATLAB/HDF5Plugin' );
    os.system( 'export HDF5_PLUGIN_PATH=~/Documents/MATLAB/HDF5Plugin' );
    
    # os.environ[ 'LD_LIBRARY_PATH' ] = "my_path" # visible in this process + all children
    
    #==================================================================================================
    
    # detect the current working directory and print it
    path = os.getcwd( );
    print( "The current working directory is %s" % path );
    
    #==================================================================================================
    
    results_path = '/home/ash/Desktop/test0_24Jan2021/results/';
    
    if not os.path.exists( results_path ):
        os.makedirs( results_path );
    
    #==================================================================================================
    
    #subfolders = ("a", "b", "c", str( 10 ) );
    subfolders = map( str, range( 0, 10 ));
    subfolders = [ s + '-run' for s in subfolders ];
    
    for subfolder in subfolders:
        #newpath = os.path.join( results_path, subfolder );
        #os.makedirs( newpath, exist_ok = True );
        os.makedirs( os.path.join( results_path, subfolder ), exist_ok = True );
   


    print( os.path.join( results_path, subfolders[ 4 ] ) )
    
    #a = ( 'cp','./*.jpg','folder_B/' ); 
    a = 'cp ' + './*.jpg ' + os.path.join( results_path, subfolders[ 4 ] + '/' ); 

    print( a )

    

    #sys.exit( 0 );



    # sys.exit( 0 );

    #for ii in range( 0, 10 ):
    #    os.replace( '*.jpg', subfolders[ ii ] );
    #    shutil.copyfile( '*.mat', dst ); 


    #==================================================================================================
    
    for ii in range( 0, 10 ):
        os.system( 'echo "zjiang202011_L0090_D15v2_006deg_run; quit" | matlab -nodisplay' );
        #os.replace( '*.jpg', subfolders[ ii ] );
        os.system('mv ' + './*.jpg ' + os.path.join( results_path, subfolders[ ii ] ) + '/' );

    #==================================================================================================
    
    
    sys.exit( 0 );



#print ("Always executed")
# 
#if __name__ == "__main__": 
#    print ("Executed when invoked directly")
#else: 
#    print ("Executed when imported")





## Python program to execute 
## function directly 
#def my_function(): 
#    print ("I am inside function")
# 
## We can test function by calling it. 
#my_function() 
#
## Python program to use 
## main for function call.
#if __name__ == "__main__":
#    my_function()
# 
#import myscript
# 
#myscript.my_function()




if __name__ == '__main__':
    try:
        main( )
    except KeyboardInterrupt:
        print( 'Interrupted' )
        try:
            sys.exit( 0 )
        except SystemExit:
            os._exit( 0 )








