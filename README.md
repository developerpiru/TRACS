# TRACS
Toolset for Ranked Analysis of CRISPR Screens - a GUI tool to analyze CRISPR screens

# Introduction
TRACS is a GUI (graphic user interface) based tool to analyze CRISPR screens. TRACS uses a ranking algorithm to identify sgRNAs and their respective genes that dropout or become enriched in experimental conditions. It requires you to provide sequencing data for a negative control conditon (cells tha do not express Cas9) and from the initial library preparation (plasmid preparation).

# Installation methods
There are several ways to install TRACS and each one is covered below. TRACS is written in Python 3.6 so will run on any operating system with Python 3.6+ installed, incluing Windows, Linux, and Mac OS. However, it relies on several dependencies which are not currently available on Windows. While this means (for the time being) you cannnot run TRACS natively on Windows, you can still use Docker (or other VM solution).

Please see the instructions below for your operating system for the best way to get started.

## Local machine or cloud\remote server (headless servers)
First you need to decide if you will run TRACS on a **local machine** with access to its GUI desktop, or if you will use a **headless\remote\cloud server**. If you plan on running TRACS on a remote\cloud server running Linux, you will first need to setup its GUI desktop in order to use TRACS. This is a one time setup process. You can also use a remote Windows or Mac OS cloud server.


As with any high-throughput sequencing analysis method, the more powerful the computer, the faster TRACS tasks will complete. We therefore recommend the use of cloud servers whenever possible, such as [Amazon AWS](https://aws.amazon.com/), [Google Cloud Platform (GCP)](https://cloud.google.com/), or [Microsoft Azure](https://azure.microsoft.com/). You may also use [Mac in Cloud](https://www.macincloud.com/), [XCloud](https://xcloud.me/), or [MacStadium](https://www.macstadium.com/) for Mac OS cloud servers.

A Linux server running Ubuntu 18.04 LTS is recommended for AWS, GCP, Azure, or other Linux VM cloud provider. If you are unfamiliar with Linux, you can also use a Windows 10 Pro instance on Azure and install TRACS using Docker. In Azure, select a [Dv3 or Ev3 based VM which support Nested Virtualization](https://azure.microsoft.com/en-us/blog/introducing-the-new-dv3-and-ev3-vm-sizes/) so that you can use Docker. At this time, AWS does not support Nested Virtualization while [GCP supports it on Linux only](https://cloud.google.com/compute/docs/instances/enable-nested-virtualization-vm-instances).

**Note:** Nested Virtualization is only required to run Linux-based Docker containers in a Linux VM. Hence if you want to do this, use Azure.


## Using TRACS natively or in a Docker container
Secondly, you need to decide whether you will run **TRACS natively** on your local or remote server, or if you will run it inside a **Docker container** (a VM instance). Our Docker container provides a ready-to-go virtual instance with TRACS and all of its required components installed. All you have to do is download and build the container (a one time step). Of course, you must have Docker isntalled. This will allow you to run TRACS on Windows, Linux, or Mac OS.


You can also install TRACS natively without the need for Docker in Linux and Mac OS. To do this, you must have Python 3.6+ installed, along with Bowtie2, Cutadapt, and MAGeCK 0.5.5.

---
## Select your desired installation method
+ Install on headless\remote\cloud server
	+ [Install Remote Desktop or VNC server & client](https://github.com/developerpiru/TRACS#install-on-headlessremotecloud-server)
	+ [Install TRACS natively](https://github.com/developerpiru/TRACS#natively-installing-tracs-currently-only-linux-and-mac-os-are-supported-natively) or [install TRACS Docker container](https://github.com/developerpiru/TRACS#using-tracs-with-docker)
+ Install on local computer
	+ [Install TRACS natively](https://github.com/developerpiru/TRACS#natively-installing-tracs-currently-only-linux-and-mac-os-are-supported-natively) or [Install TRACS Docker container](https://github.com/developerpiru/TRACS#using-tracs-with-docker)
---

# Install on headless\remote\cloud server

You must first connect to your server's GUI desktop. If you are using a Windows server (Windows 10 or any Windows Server OS), you can connect to it using Remote Desktop Connection from any host Windows OS and you will see the desktop automatically. If you are using a Mac OS server, you can you can install the Microsoft Remote Desktop App through the App Store to connect to a remote Mac OS or Windows server. Alternatively, you may follow the instructions below to setup VNC on your Mac OS server. If you are connecting to a Linux server, you must first ssh into your server and setup a VNC server as described below.

If you are using a Linux server, chances are it is currently only configured for ssh\terminal access. Since TRACS is a GUI program, it is critical that you have an X window system setup on your Linux server. **If you already have this setup and you can connect to the Ubuntu desktop using VNC, you can skip this part and proceed to the instructions for installing TRACS natively.**

## Setting up VNC on Ubuntu (may also work for Mac OS but untested)
1. Install VNC server and xfce4 components:
	```
	sudo apt update
	sudo apt install xfce4 xfce4-goodies
	```
2. Install Tight VNC Server:
	```
	sudo apt install tightvncserver
	```
3. Run initial VNC server configuration to setup a password. You can say no to a "view only password":
	```
	vncserver
	```
4. Exit VNC server:
	```
	vncserver -kill :1
	```
5. Edit the VNC server configuration file to start xfce4 desktop on launch:
	```
	sudo nano ~/.vnc/xstartup
	```
6. Replace the contents of the file with the following:
	```
	#!/bin/bash
	xrdb $HOME/.Xresources
	startxfce4 &
	```
7. Save the file and exit (press Control + X and confirm to save).

8. Edit the permissions of the file:
	```
	sudo chmod +x ~/.vnc/xstartup
	```
	
## Install VNC client on your local computer
You need a VNC client to connect to the VNC server you just setup on your remote server. 

You can install either PuTTY for Windows (https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html) or RealVNC's VNC Viewer for Windows, Mac OS, or Linux (https://www.realvnc.com/en/connect/download/viewer/)

## Connecting to VNC server on remote server
1. Start the VNC server with this command on your remote server:
	```
	vncserver -geometry 1200x1050
	```
	**Note:** the `-geometry 1200x1050` flag is optional and you can customize it to any resolution you prefer

2. Open your VNC client (PuTTY or VNC Viewer) and enter the **external IP address** (not its local IP!) of your remote computer followed by the `:5901` port. For example: `192.0.2.0:5901`
	
3. Enter the password you created earlier when prompted.

You should now see the Ubuntu desktop now and be able to interact with the GUI using your mouse and keyboard. 

You are now connected to your remote server's desktop interface and can continue with the installation for TRACS either [natively](https://github.com/developerpiru/TRACS#natively-installing-tracs-currently-only-linux-and-mac-os-are-supported-natively) or [using Docker](https://github.com/developerpiru/TRACS#using-tracs-with-docker).

**Note:** to help you with downloading required components, it is recommended that you install Chrome or Firefox on your Linux server:

Chrome:
	
	wget https://dl.google.com/linux/direct/google-chrome-stable_current_amd64.deb
	sudo dpkg -i google-chrome-stable_current_amd64.deb
	
	
Firefox:
	
	sudo apt-get update
	sudo apt-get install firefox
	

#### Stopping the VNC server
You can stop the VNC server using this command:
	
	vncserver -kill :1
	
---
## Natively installing TRACS (currently only Linux and Mac OS are supported natively)

Once you can connect to your remote server's GUI desktop, or if you are installing TRACS on a local machine, you can proceed with the following instructions to install TRACS natively. If you'd like to install TRACS in a Docker container instead, [please see these instructions](https://github.com/developerpiru/TRACS#using-tracs-with-docker).

### Requirements
You must have the following components installed on your Linux or Mac OS device to run TRACS:
	1. Python 3.6+
	2. Tkinter, pandas, and numpy packages for Python 3
	3. Cutadapt (recommended version 1.18)
	4. Bowtie2 (recommended version 2.3.4.2)
	5. Samtools (recommended version 1.9)
	6. MAGeCK (only version 0.5.5 supported)

If you try to run TRACS and get errors during launch, it is likely you are missing the Tkinter, Pandas, or Numpy packages.

If you get errors during analysis, you likely don't have Cutadapt, Bowtie2, Samtools, or MAGeCK properly installed or configured.

**Note:** If you are using Ubuntu 18.04 LTS or later, chances are you already have Python 3.6+ installed.

If you are unsure whether you have these components installed, here is an easy way to check (run these commands in a Terminal:

**Python 3.6+:**
	
	python --version
		
If installed, you should see the current version of Python.


**Cutadapt:**

	cutadapt --version
	
If installed, you should see the current version of Cutadapt.


**Bowtie2:**

	bowtie2 --version
	
If installed, you should see the current version of Bowtie2.


**Samtools:**

	samtools --version
	
If installed, you should see the current version of Samtools.


**MAGeCK:**

	mageck --version
	
If installed, you should see the current version of MAGeCK.
  
  
To determine if you have the correct **Python packages installed**, start Python:
	
	python3
	
Then import each package from the Python 3 command prompt:
	
	>>>  import tkinter

	>>>  import pandas
	
	>>>  import numpy
	
	>>>  exit()


If you are able to import each without errors, you have them installed. If you receive a `ModuleNotFoundError` error, you need to install that package.

If you have each of these required components installed, you can continue with the instructions to [natively installing TRACS in Linux or Mac OS](https://github.com/developerpiru/TRACS#installing-tracs-natively).

If these required components are missing or you don't have at least the recommended versions installed, you can follow the [instructions below to install them](https://github.com/developerpiru/TRACS#installing-required-components). 

**Note:** For MAGeCK, we find versions after 0.5.5 cause errors. For the purposes of TRACS, MAGeCK is only used to generate read counts so version 0.5.5 is sufficient and need not be upgraded.

### Installing required components
Update sources before you begin:
	
	sudo apt-get update
	
	
Installing Python 3:

	sudo apt-get install python3


Installing tkinter for Python3:

	sudo apt-get install python3-tk

	
Installing pandas for Python3:

	sudo apt-get install python3-pandas

Installing pandas will usually install numpy as well. In case it didn't, you can install it as follows:

	sudo apt-get install python3-numpy


Installing Cutadapt:

	sudo apt-get install cutadapt


Installing Bowtie2:

	sudo apt-get install cutadapt


Installing Samtools:

	sudo apt-get install samtools 

	
Installing MAGeCK:
To install MAGeCK, we recommend building it from source. This requires some additional packages not required by TRACS but required to build from source:

	sudo apt-get install python3-distutils python3-dev build-essential


Then download MAGeCK 0.5.5 using ```wget``` or paste the link into a browser:

	wget http://www.dkfz.de/signaling/crispranalyzer/mageck-0.5.5.tar.gz

	
Then extract and run the setup script (make sure you are in the correct directory if you used a browser to download!):

	tar xvf mageck-0.5.5.tar.gz
	cd mageck-0.5.5
	python setup.py install


Go back to the previous steps under the [Requirements](https://github.com/developerpiru/TRACS#requirements) section and ensure each component is installed as described above. Otherwsie, you can now [install TRACS natively](https://github.com/developerpiru/TRACS#installing-tracs-natively).

### Installing TRACS natively
Once you have all of the requirement components installed, you are ready to install TRACS!

1. Access the desktop of your Linux or Mac OS device. You need to connect remotely to the desktop if you are using a remote\cloud server ([see instructions above](https://github.com/developerpiru/TRACS#install-vnc-client-on-your-local-computer)).

2. Download the TRACS zip file from our github repository ([TRACSv0.9.9.2.zip](TRACSv0.9.9.2.zip)).

3. Extract the TRACSv0.9.9.2 folder to any location you desire.

4. Open a Terminal window and nagvigate to that folder. 

5. Enter this command to start TRACS:
	```
	python3 TRACS_v0.9.9.2.py
	```
---
## Using TRACS with Docker

Access the desktop of your Windows, Linux or Mac OS device. You need to connect remotely to the desktop if you are using a remote\cloud server ([see instructions above](https://github.com/developerpiru/TRACS#install-on-headlessremotecloud-server)).

### Docker on Windows
#### Installation
1. Download and install Docker Desktop for Windows 10 Pro\Enterprise here: https://www.docker.com/products/docker-desktop

	If you have an older version of Windows (or Windows 10 Home) you will need to install Docker Toolbox: https://docs.docker.com/toolbox/toolbox_install_windows/

	**Note:** you may have to sign up for free with Docker before you can download.

	Follow the instructions on the respective Docker pages for installing Docker and getting it running.

2. Download and install XMing X Server for Windows: https://sourceforge.net/projects/xming/
	
	Install with the default options.

3. Download the [TRACS-docker-container-win.zip](TRACS-Docker-setup/TRACS-docker-container-win.zip) file located in the TRACS-Docker-Setup folder on our github repository. 
	
	This contains the TRACS Dockerfile, XMing configuration, and Windows PowerShell scripts to automate setup and launching TRACS. 

4. Extract the TRACS-docker-container folder to ```C:\``` (note if you change this location, you will need to modify the PowerShell scripts accordingly).

5. Right click the "Build TRACS Container.ps1" file and click "Run with PowerShell". 
	Wait for the process to complete; it may take several minutes. 
	
#### Launching TRACS Docker container on Windows
Right click the "Run TRACS Container - with XLaunch.ps1" file and click "Run with PowerShell". 
Approve the Windows access control prompt if necessary and approve the sharing of your local ```C:\``` if prompted by Docker.

TRACS will launch in a Docker container and mount your local ```C:\``` drive at ```/app/TRACS/cdrive/``` in the Docker container so you can transport files from the container to your local drive. 

**IMPORTANT!:** As with any Docker container, anything you DO NOT save in ```/app/TRACS/cdrive/``` will be lost when you exit TRACS!

---
### Docker on Mac OS
#### Installation
1. Download and install Docker Desktop for Mac OS here: https://www.docker.com/products/docker-desktop

	If you have an older version of Mac OS  you will need to install Docker Toolbox: https://docs.docker.com/toolbox/toolbox_install_mac/

	**Note:** you may have to sign up for free with Docker before you can download.

	Follow the instructions on the respective Docker pages for installing Docker and getting it running.

2. Download and install XQuartz for Mac OS: https://www.xquartz.org/
	
	Install with the default options.

3. Start XQuartz and go to Preferences > Security and check option to allow connections from network clients. 
	
	Exit XQuartz.

4. Download the [TRACS-docker-container-mac.zip](TRACS-Docker-setup/TRACS-docker-container-mac.zip) file located in the TRACS-Docker-Setup folder on our github repository. 
	
	This contains the TRACS Dockerfile, XMing configuration, and Windows PowerShell scripts to automate setup and launching TRACS. 

5. Extract the TRACS-docker-container folder to anywhere you desire.

6. Double click the "Build-TRACS" file to start building the TRACS Docker container.

#### Launching TRACS Docker container on Mac OS
Double click the "Start-TRACS" file to start the Docker container and launch TRACS. 

TRACS will launch in a Docker container and mount your local drives (```/Volumes```) at ```/app/TRACS/LocalDrives/``` in the Docker container so you can transport files from the container to your local drive. 

**IMPORTANT!:** As with any Docker container, anything you DO NOT save in ```/app/TRACS/LocalDrives/``` will be lost when you exit TRACS!

---
### Docker on Linux
#### Installation
1. Install Docker for Linux (Ubuntu):
	```
	sudo apt-get update
	sudo apt-get install docker.io
	```
	
2. Start Docker:
	```
	sudo systemctl start docker
	sudo systemctl enable docker
	
	```
3. Download the [TRACS-docker-container-linux.zip](TRACS-Docker-setup/TRACS-docker-container-linux.zip) file located in the TRACS-Docker-Setup folder on our github repository. 

4. Extract the TRACS-docker-container folder to anywhere you desire:
	```
	unzip TRACS-docker-container-linux.zip
	```

5. Navigate to the folder and build the TRACS Docker container:
	```
	cd TRACS-docker-container
	docker build -t tracs .
	```
	
#### Launching TRACS Docker container on Linux
Open a terminal window and enter this command to start the TRACS container:
```
docker run -ti --rm -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix ubuntu -v /path/to/folder:/app/TRACS/sharedfolder tracs
```
			
Where ```/path/to/folder``` is your local drive/folder that you want to make available to the TRACS container.
	
For example:
```
docker run -ti --rm -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix ubuntu -v $HOME/:/app/TRACS/sharedfolder tracs
```

TRACS will launch in a Docker container and mount your local drive or folder (```/path/to/folder```) at ```/app/TRACS/sharedfolder/``` in the Docker container so you can transport files from the container to your local drive.

**IMPORTANT!:** As with any Docker container, anything you DO NOT save in ```/app/TRACS/sharedfolder/``` will be lost when you exit TRACS!

---

# Demo files

Check the [Demo-files.md](Demo-files.md) file for information on how to download some demo files for testing TRACS.
