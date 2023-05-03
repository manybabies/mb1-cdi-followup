# Reproducible Research
The ability to reproduce results of computational research is becoming a critical issue with the onging technicial changes and computational updates. Docker, an open source technology, addresses this issue  by providing a virtual environment with its package and dependencies to ensure a long-term reproducibilty of scientific results. In this tutorial, we will walk you through a few steps to install Docker software on your local machine, and to import and run *mb1-cdi-followup*'s Rstudio Docker image.  A *Docker image* is a pre-configured copy of software, Rstudio in our case, with an installation of selected R versions, packages, and dependencies. In other words, a *Docker image* is a virtual copy of *mb1-cdi-followup*'s original R environment that was used to produce the results. 

In order to run *mb1-cdi-followup*'s Rstudio Docker image, you have to follow three steps.
- Step 1: Create Docker Account.
- Step 2: Download and Install Docker Software.
- Step 3: Import and Run *mb1-cdi-followup*'s Rstudio Docker Image.


## Step 1: Create Docker Account
- Create an account in [Docker](https://hub.docker.com/). Please, keep a copy of your *account ID* and *password* to use them in *step 3*.

## Step 2: Download & Install Docker Software
- [Docker for Mac](#docker-for-mac)
- [Docker for W

ows](#docker-for-windows)


### Docker for Mac
**Note**: Keep in mind that your macOS must be at least 10.15 or newer (*Catalina*, *Big Sure*, *Monterey*)  to be able to run Docker.

#### Insallation: 
1. [Docker for Mac](https://desktop.docker.com/mac/main/amd64/Docker.dmg?utm_source=docker&utm_medium=webreferral&utm_campaign=docs-driven-download-mac-amd64).
2. Double-click the `Docker.dmg` file, then drag-and-drop `Docker.app` in your Application folder.
3. Open `Docker.app` to start the application, and sign in  with your Docker account ID and passcode.
4. Locate the "whale" icon on your status bar, then click it. Be sure that it indicates `Docker is running` as shown below.

<p align="center">
  <img width="200" src="https://user-images.githubusercontent.com/47132064/184414404-529568c9-f577-4787-a372-e8e14431b896.png">
</p>
5. Some macOS such as *Catalina*  might have an issue with Docker by showing the persistent message as `Docker is starting` on the whale icon. In this case, please open Docker application and click on troublshouing icon as shown below. Then, click on `Reset to factory default`. Next, you have to restart the application. It should work fine afterwards.

<p align="center">
  <img width="600" src="https://user-images.githubusercontent.com/47132064/184058430-5f7c6af4-deeb-4d92-9be3-07efd9ea9174.jpg">
</p>

- If the whale icon indicates `Docker is running`, then proceed to **Step 3**.

### Docker for Windows
**Note**: Windows machines must have a 64-bit processor and at least 4GB system RAM to successfully run Docker on Windows. Follow the guidelines in this [link](https://support.microsoft.com/en-gb/topic/determine-whether-your-computer-is-running-a-32-bit-version-or-64-bit-version-of-the-windows-operating-system-1b03ca69-ac5e-4b04-827b-c0c47145944b) to determine if your machine supports 64-bit processors.
#### Insallation: 

1. [Download Docker for Windows](https://desktop.docker.com/win/main/amd64/Docker%20Desktop%20Installer.exe), then double-click on the `Docker Desktop Insaller.exe` file to start the installation process.
2. Once the installation is successful, open the Docker software. After a few seconds, it should indicate `docker is running`.  Then, proceed to *step 3*.



## Step 3: Import and run *mb1-cdi-followup*'s Rstudio Docker Image
From your local machine Open the following application:
  - For **MAC**, open `Terminal`.
  - For **Windows**, open `Command Prompt`.
  
Then, follow the instructions below. For each time, copy the command line without the dollar sign `$`, and paste it in `Terminal` or `Command Prompt`.
1. Log in with your Docker account ID and password: `$ docker login`
  - Enter your Docker account ID.
  - Enter your password. **Note**:  while you are typing your password, it will not show up and the cursor will not move. Just type it, then hit *enter*/*return*.
  
  ( In your next log in, you may not your need to enter account ID and password. Once you enter the login command line, it will log you in automatically.)
  
2. Then, import *mb1-cdi-followup*'s Rstudio Docker image: `$ docker pull mb1l/mb1-cdi-followup:v.1.2`
3. Next, run *mb1-cdi-followup*'s Rstudio Docker image to be able to reproduce the results: `$ docker run -p 8787:8787 -e PASSWORD=rstudio mb1l/mb1-cdi-followup:v.1.2`
4. Once the `Terminal`/`Command Prompt` indicates `[services.d] done.`, open your browser, then copy & paste the following: `localhost:8787`.

- Copy & paste the following *user ID* & *password*:
  - User ID: `rstudio`.
  - Password: `rstudio`.
- Once you lgo in, Rstudio should look like the image below:
<p align="center">
  <img width="800" src="https://user-images.githubusercontent.com/47132064/184929809-4888c209-e074-4446-be25-e06b1ce75487.jpg">
</p>


## Stop Docker Image

After finishing reproducing the results, you should stop the running image. In  `Terminal` or `Command Prompt` pressing `Control (CTRL) + c`. Afterward, you should get the following message `[finish] sending all processes the KILL signal and exiting`.

## Remove Docker

You can reclaim the used storage in your machine after reproducing the results. Follow these steps:
1. In  `Terminal` or `Command Prompt` enter the following command line: `$ docker system prune -a`.
2. Uninstall Dokcer desktop from your `application folder` or `program files`

