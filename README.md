# TRACS
**T**oolset for the **R**anked **A**nalysis of **C**RISPR **S**creens - a GUI tool to analyze CRISPR screens

# Introduction
TRACS is a GUI (graphic user interface) based tool to analyze CRISPR screens. TRACS uses a ranking algorithm to identify sgRNAs and their respective genes that dropout or become enriched in experimental conditions. It requires you to provide sequencing data for a negative control conditon (cells tha do not express Cas9) and from the initial library preparation (plasmid preparation). 

Data output from TRACS can be visualized using the companion app, [VisualizeTRACS](https://github.com/developerpiru/VisualizeTRACS).

---

# Table of contents
+ [Using TRACS with Docker](https://github.com/developerpiru/TRACS#using-tracs-with-docker)
	+ [Docker on Windows](https://github.com/developerpiru/TRACS#docker-on-windows)
	+ [Docker on Mac OS](https://github.com/developerpiru/TRACS#docker-on-mac-os)
	+ [Docker on Linux](https://github.com/developerpiru/TRACS#docker-on-linux)
+ [Demo files](https://github.com/developerpiru/TRACS#demo-files)
+ [Data visualization and exploration](https://github.com/developerpiru/TRACS#step-5-data-visualization-and-exploration)
---


## Using TRACS with Docker

Connect to the desktop of your Windows, Linux or Mac OS device. You need to connect remotely to the desktop if you are using a remote\cloud server ([see instructions above](https://github.com/developerpiru/TRACS#install-on-headlessremotecloud-server)).

- [Docker on Windows](https://github.com/developerpiru/TRACS#docker-on-windows)
- [Docker on Mac OS](https://github.com/developerpiru/TRACS#docker-on-mac-os)
- [Docker on Linux](https://github.com/developerpiru/TRACS#docker-on-linux)

### Docker on Windows
#### Installation
1. Download and install Docker Desktop for Windows 10 Pro\Enterprise here: https://www.docker.com/products/docker-desktop

	If you have an older version of Windows (or Windows 10 Home) you will need to install Docker Toolbox: https://docs.docker.com/toolbox/toolbox_install_windows/

	Follow the instructions on the respective Docker pages for installing Docker and getting it running.

2. Download and install XMing X Server for Windows: https://sourceforge.net/projects/xming/
	
	Install with the default options.

3. Download the [TRACS-docker-launcher-windows.zip](TRACS-Docker-setup/TRACS-docker-launcher-windows.zip) file located in the TRACS-Docker-Setup folder on our github repository. 
	
	This contains the launcher for TRACS, XMing configuration, and Windows PowerShell scripts to automate setup and launching TRACS. 

4. Extract the TRACS-docker-launcher-windows folder. 

5. Right click the "Launch TRACS.ps1" file and click "Run with PowerShell". 
	Wait for the process to complete; it may take several minutes.
	The script will pull the latest TRACS image from Docker, use the XMing configuration script, and launch the TRACS interface.
	Approve the Windows access control prompt if necessary and approve the sharing of your local ```C:\``` if prompted by Docker.

TRACS will launch from the Docker container and mount your local ```C:\``` drive at ```/app/TRACS/cdrive/``` in the Docker container so you can transport files from the container to your local drive. 

**IMPORTANT!:** As with any Docker container, anything you DO NOT save in ```/app/TRACS/cdrive/``` will be lost when you exit TRACS!

---
### Docker on Mac OS
#### Installation
1. Download and install Docker Desktop for Mac OS here: https://www.docker.com/products/docker-desktop

	If you have an older version of Mac OS  you will need to install Docker Toolbox: https://docs.docker.com/toolbox/toolbox_install_mac/

	Follow the instructions on the respective Docker pages for installing Docker and getting it running.

2. Download and install XQuartz for Mac OS: https://www.xquartz.org/
	
	Install with the default options.

3. Start XQuartz and go to Preferences > Security and check option to allow connections from network clients. 
	
	Exit XQuartz.

4. Download the [Start-TRACS-MacOS](TRACS-Docker-setup/Start-TRACS-MacOS) file located in the TRACS-Docker-Setup folder on our github repository. 

5. Double click the "Start-TRACS" file to start the Docker container and launch TRACS.
	Wait for the process to complete; it may take several minutes.
	The script will pull the latest TRACS image from Docker, use the XMing configuration script, and launch the TRACS interface.

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

5. Open a terminal window and enter this command to start the TRACS container:

	docker run -ti --rm -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix ubuntu -v /path/to/folder:/app/TRACS/sharedfolder pirunthan/tracs:latest

			
Where ```/path/to/folder``` is your local drive/folder that you want to make available to the TRACS container.
	
For example:
	
	docker run -ti --rm -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix ubuntu -v $HOME/:/app/TRACS/sharedfolder pirunthan/tracs:latest
	

TRACS will launch in a Docker container and mount your local drive or folder (```/path/to/folder```) at ```/app/TRACS/sharedfolder/``` in the Docker container so you can transport files from the container to your local drive.

**IMPORTANT!:** As with any Docker container, anything you DO NOT save in ```/app/TRACS/sharedfolder/``` will be lost when you exit TRACS!

---

# Demo files

Check the [Demo-files.md](Demo-files.md) file for information on how to download some demo files for testing TRACS.

---

# Data visualization and exploration

You can visualize the data in your ```[ExperimentName].csv``` TRACS output file using the companion app called [VisualizeTRACS](https://github.com/developerpiru/VisualizeTRACS), which is an R shiny app. Unlike TRACS, you do not need to have a powerful computer to run VisualizeTRACS - this means you can download the ```[ExperimentName].csv``` file from a cloud server (if you used one) and do the visualization steps locally. 
