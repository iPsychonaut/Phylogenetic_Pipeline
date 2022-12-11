import platform
import subprocess
import urllib.request
import os
import zipfile
import sys
import subprocess

# URL of the repository to download
PhyloPipeline_url = 'https://github.com/iPsychonaut/Phylogenetic_Pipeline/archive/refs/heads/main.zip'

# Get the path to the user's home directory
home_path = os.path.expanduser('~')

# Get the path to the user's main drive
base_path = home_path.split('Users')[0] 

# Get the path to the user's downloads folder
downloads_path = f'{home_path}/Downloads'

# Destination path where the program will be installed
install_path = f'{base_path}PhyloPipeline'

# Destination path where the dependencies and images will be downloaded to
bin_path = f'{install_path}/bin'

# Use the platform module to get the system name and version
system_name = platform.system()
system_version = platform.version()

# Print the system name and version
print(f'Operating system: {system_name}')

# Build out link lists for each EXE
windows_exe_links = ['https://github.com/Cibiv/IQ-TREE/releases/download/v1.6.12/iqtree-1.6.12-Windows.zip',
                     'http://trimal.cgenomics.org/_media/trimal.v1.2rev59.zip',
                     'https://github.com/rcedgar/muscle/releases/download/5.1.0/muscle5.1.win64.exe',
                     'https://github.com/rambaut/figtree/releases/download/v1.4.4/FigTree.v1.4.4.zip']

mac_exe_links = ['https://github.com/Cibiv/IQ-TREE/releases/download/v1.6.12/iqtree-1.6.12-MacOSX.zip',
                 'http://trimal.cgenomics.org/_media/trimal.v1.2rev59.tar.gz']

linux_exe_links = ['https://github.com/Cibiv/IQ-TREE/releases/download/v1.6.12/iqtree-1.6.12-Linux.tar.gz',
                   'http://trimal.cgenomics.org/_media/trimal.v1.2rev59.tar.gz']


# Download the repository as a zip file
PhyloPipeline_zip_path = f'{downloads_path}/repository.zip'
urllib.request.urlretrieve(PhyloPipeline_url, PhyloPipeline_zip_path)

# Extract the zip file
with zipfile.ZipFile(PhyloPipeline_zip_path, 'r') as zip:
    zip.extractall(downloads_path)

# Move the extracted folder to the destination folder with the desired name
os.rename(f'{downloads_path}/Phylogenetic_Pipeline-main', install_path)

# Remove Downloaded Files
os.remove(PhyloPipeline_zip_path)

def figtree_install():
    if system_name == 'Windows':
        # Check if EXE is in the bin directory
        isExist = os.path.exists(f'{bin_path}/FigTree v1.4.4.exe')
        if not isExist:
            figtree_url = exe_links[3]
            
            # Set the name of the downloaded file
            figtree_zip = f'{downloads_path}/{figtree_url.split("/")[-1]}'
            figtree_zip_name = figtree_zip.split("/")[-1].replace(".zip","").replace('.v',' v')
            # Download the file from the URL
            urllib.request.urlretrieve(figtree_url, figtree_zip)
            
            #Extract the contents of the zip file to the current directory
            with zipfile.ZipFile(figtree_zip, 'r') as zip:
                zip.extractall(bin_path)
        
            os.rename(f'{bin_path}/{figtree_zip_name}',f'{bin_path}/figtree')
            os.rename(f'{bin_path}/figtree/{figtree_zip_name}.exe',f'{bin_path}/figtree/figtree.exe')
            # Remove Downloaded Files
            os.remove(figtree_zip)
            print(f'FigTree Successfully Installed!')
        else:
            print('FigTree already installed.')
    elif system_name == 'Darwin':
        print('DO MAC FIGTREE')
    elif system_name == 'Linux':
            print('DO LINUX FIGTREE')

def iqtree_install():

    # Check if IQ-Tree is in the bin directory
    isExist = os.path.exists(f'{bin_path}/iqtree.exe')
    if not isExist:
        iqtree_url = exe_links[0]
        
        # Set the name of the downloaded file
        iqtree_zip = f'{downloads_path}/{iqtree_url.split("/")[-1]}'
       
        # Download the file from the URL
        urllib.request.urlretrieve(iqtree_url, iqtree_zip)
        
        #Extract the contents of the zip file to the current directory
        with zipfile.ZipFile(iqtree_zip, 'r') as zip:
            zip.extractall(bin_path)
        
        # Move the EXE file to the bin folder
        os.rename(f'{bin_path}/{iqtree_url.split("/")[-1].replace(".zip","")}', f'{bin_path}/iqtree')  
        
        # Remove Downloaded Files
        os.remove(iqtree_zip)
        print('IQ-Tree Successfully Installed!')
    else:
        pass
        print('IQ-Tree is already installed.')

def trimAl_install():

    # Check if TrimAl is in the bin directory
    isExist = os.path.exists(f'{bin_path}/trimal.exe')
    if not isExist:
        trimAl_url = exe_links[1]
        
        # Set the name of the downloaded file
        trimAl_zip = f'{downloads_path}/{trimAl_url.split("/")[-1]}'
       
        # Download the file from the URL
        urllib.request.urlretrieve(trimAl_url, trimAl_zip)
        
        # Open the zip file
        with zipfile.ZipFile(trimAl_zip, 'r') as zip:
            
            # Extract the contents of the zip file to the current directory
            zip.extractall(bin_path)
        
        # Remove Downloaded Files
        os.remove(trimAl_zip)
        print('TrimAl Successuflly Installed!')
    else:
        pass
        print('TrimAl is already installed.')

def conda_check():
    try:
        import conda # Try importing the 'conda' module
    
        # If the most recent version of Anaconda is not installed, install it using the 'conda' command
        subprocess.run(['conda', 'install', 'anaconda'])
    
        # Print a message to confirm that Anaconda was installed
        print('The most recent version of Anaconda has been installed.')
    except ImportError:
        
        # If Anaconda is not installed, install it using the 'conda' command
        try:
            subprocess.run(['conda', 'install', 'anaconda'])
        
            # Print a message to confirm that Anaconda was installed
            print('Anaconda Successfully Installed!')
        except FileNotFoundError:
            print('Anaconda already installed.')

def muscle_install():
    if system_name == 'Windows':
        
        # Set the URL to download from
        muscle_url = exe_links[2]
        
        # Set the name of the downloaded file
        muscle_exe = f'{bin_path}/muscle.exe'
         
        # Download the file from the URL
        urllib.request.urlretrieve(muscle_url, muscle_exe)
    else:
        # If the library is not installed, install it using the 'conda' command
        subprocess.run(['conda', 'install', 'muscle'])

    # Print a message to confirm that the library was installed
    print('MUSCLE Successfully installed!')


def python_check():
    # Check if the 'python' command is available on the system
    if 'python' in sys.executable:
        # If the 'python' command is available, print a message
        print('Python is already installed.')
    else:
        # If the 'python' command is not available, Use a 'case' statement to
        # select the appropriate package manager or installer based on the operating system
        match system_name:
    
            # For Linux systems, use the 'apt-get' command
            case 'Linux':
                subprocess.run(['apt-get', 'install', 'python3'])
    
            # For Windows systems, use the official Python installer
            case 'Windows':
                
                # Download from the offical website
                python_url = 'https://www.python.org/ftp/python/3.11.1/python-3.11.1-amd64.exe'
    
                # Set the name of the downloaded file
                python_exe = python_url.split('/')[-1]
         
                # Download the file from the URL
                urllib.request.urlretrieve(muscle_url, python_exe)
                subprocess.run([python_exe])
    
            # For Mac systems, use the 'brew' command
            case "Darwin":
                subprocess.run(['brew', 'install', 'python'])
    
            # For other operating systems, print a message
            case _:
                print('Unsupported operating system.')
    
        # Print a message to confirm that Python was installed
        print('Python Successfully Installed!')
    
    # Check whether the specified path exists or not
    isExist = os.path.exists(bin_path)
    if not isExist:
        
        # Create a new directory because it does not exist
        os.mkdir(bin_path)

# RUN MAIN CODE
python_check()
conda_check()
if system_name == 'Windows':
    exe_links = windows_exe_links
    muscle_install()
    trimAl_install()
    iqtree_install()
    figtree_install()
elif system_name == 'Darwin':
    exe_links = mac_exe_links
    muscle_install()
    trimAl_install()
    iqtree_install()
    figtree_install()
elif system_name == 'Linux':
    exe_links = linux_exe_links
    muscle_install()
    trimAl_install()
    iqtree_install()
    figtree_install()
else:
    print(f'NO DEPENDANCY DOWNLOADS AVAILABLE FOR {system_name}')