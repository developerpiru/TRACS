# TRACS
Toolset for Ranked Analysis of CRISPR Screens - a GUI tool to analyze CRISPR screens

# Introduction
TRACS is a GUI (graphic user interface) based tool to analyze CRISPR screens. TRACS uses a ranking algorithm to identify sgRNAs and their respective genes that dropout or become enriched in experimental conditions. It requires you to provide sequencing data for a negative control conditon (cells tha do not express Cas9) and from the initial library preparation (plasmid preparation).

# Installation
There are several ways to install TRACS are each one is covered below. TRACS is written in Python 3.6 so will run on any operating system with Python 3.6+ installed, incluing Windows, Linux, and Mac OS. However, it relies on several dependencies which are not currently available on Windows. Please see the instructions below for your operating system for the best way to get started.

## Local machine or cloud\remote server (headless servers)
First decide if you will run TRACS on a local machine with access to its graphical desktop or if you will use a remote server\computer. If you already have a computer that you can access locally and get to the operating system's GUI desktop, continue below to see instructions for that respective operating system. If you plan on running TRACS on a remote\cloud server running Linux, you will first need to setup its GUI desktop in order to use TRACS. This is a one time setup process.

## Using TRACS natively or running a Docker container
The easiet way to install TRACS (and all of its required components) is to use our Docker container. This container image contains all of the required components to run TRACS and only requires that you have Docker installed on your operating system (Windows, Linux, or Mac OS).

## Using Docker on Windows
### Installation
1. Download and install Docker Desktop for Windows 10 Pro\Enterprise here: https://www.docker.com/products/docker-desktop

	If you have an older version of Windows (or Windows 10 Home) you will need to install Docker Toolbox: https://docs.docker.com/toolbox/toolbox_install_windows/

	Note: you may have to sign up for free with Docker before you can download.

	Follow the instructions on the respective Docker pages for installing Docker and getting it running.

2. Download and install XMing X Server for Windows: https://sourceforge.net/projects/xming/
	Install with the default options.

3. Download our TRACS-docker-container-win.zip from our github repository. 
	This contains the TRACS Dockerfile, XMing configuration, and Windows PowerShell scripts to automate setup and launching TRACS. 

4. Extract the TRACS-docker-container folder to C:\ (note if you change this location, you will need to modify the PowerShell scripts accordingly).

5. Right click the "Build TRACS Container.ps1" file and click "Run with PowerShell". 
	Wait for the process to complete; it may take several minutes. 
	
### Launching TRACS Docker container on Windows
Right click the "Run TRACS Container - with XLaunch.ps1" file and click "Run with PowerShell". 
Approve the Windows access control prompt if necessary and approve the sharing of your local C:\ if prompted by Docker.

TRACS will launch in a Docker container and mount your local C:\ drive at /app/TRACS/cdrive/ in the Docker container so you can transport files from the container to your local drive. Note that as with any Docker container, anything you DO NOT save in /app/TRACS/cdrive/ will be lost when you exit TRACS!

## Using Docker on Mac OS
### Installation
1. Download and install Docker Desktop for Mac OS here: https://www.docker.com/products/docker-desktop

	If you have an older version of Mac OS  you will need to install Docker Toolbox: https://docs.docker.com/toolbox/toolbox_install_mac/

	Note: you may have to sign up for free with Docker before you can download.

	Follow the instructions on the respective Docker pages for installing Docker and getting it running.

2. Download and install XQuartz for Mac OS: https://www.xquartz.org/
	Install with the default options.

3. Start XQuartz and go to Preferences > Security and check option to allow connections from network clients. Exit XQuartz.

4. Download our TRACS-docker-container-mac.zip from our github repository. 
	This contains the TRACS Dockerfile, XMing configuration, and Windows PowerShell scripts to automate setup and launching TRACS. 

5. Extract the TRACS-docker-container folder to anywhere you desire.

6. Double click the "Build-TRACS" file to start building the TRACS Docker container.

### Launching TRACS Docker container on Mac OS
Double click the "Start-TRACS" file to start the Docker container and launch TRACS. 

TRACS will launch in a Docker container and mount your local drives (/Volumes) at /app/TRACS/LocalDrives/ in the Docker container so you can transport files from the container to your local drive. Note that as with any Docker container, anything you DO NOT save in /app/TRACS/LocalDrives/ will be lost when you exit TRACS!

## Using Docker on Linux
To be done...
