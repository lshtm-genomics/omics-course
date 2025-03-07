# Getting connected

### Introduction

There are three comonents you will need to understand before we get connected.

####  GitHub Codespaces

GitHub Codespaces is a cloud-based development environment that allows users to run and edit code directly in a virtual machine (VM) hosted by GitHub. It provides a fully configured workspace, pre-installed with necessary tools like Git, Python, and Docker, eliminating the need for users to set up a local environment.

#### Visual Studio Code (VS Code)

VS Code is a lightweight yet powerful code editor that supports multiple programming languages. It also allows users to connect to GitHub Codespaces, allowing users to work on cloud-based VMs as if they were local files.
Users will launch GitHub Codespaces via VS Code, enabling them to interact with their virtual machine through a familiar coding interface.

#### TigerVNC

TigerVNC (Virtual Network Computing) is a tool that allows users to remotely access a graphical desktop environment of their virtual machine. Since GitHub Codespaces primarily provides a terminal interface, TigerVNC will let users interact with a full virtual desktop, running graphical applications as if they were using a local computer.

#### How These Work Together:

1. Users launch GitHub Codespaces to create a cloud-based virtual machine.
2. They connect to the VM using VS Code, accessing the terminal and code files.
3. If a graphical interface is needed, they use TigerVNC to access the full desktop environment.

This setup allows users to work on cloud-hosted projects with both command-line and GUI access, making it ideal for bioinformatics and data analysis workflows.


## Getting set up

### 1. Set up a GitHub account

Head over to https://github.com/ and sign up for an account. You can skip this step if you already have an account.

### 2. Download Visual Studio Code

Get the latest version of vscode from https://code.visualstudio.com/ and follow install instructions from the website.

### 3. Download TigerVNC 

Download TigerVNC from https://sourceforge.net/projects/tigervnc/files/stable/1.15.0/ and follow install instructions from the website. Based on your operating system, you can download the appropriate version.

- Windows: tigervnc64-1.15.0.exe
- Mac: TigerVNC-1.15.0.dmg

Follow the instructions provided to install TigerVNC on your system.

### 4. Get connected

Once you have installed TigerVNC, you can connect to your GitHub Codespace using the following steps:

1. Open Visual Studio Code.
2. Open the command palette by pressing `Ctrl+Shift+P` (Windows/Linux) or `Cmd+Shift+P` (Mac).
3. Type `Create new Codespace` and select the option.
4. Enter `jodyphelan/teaching-codespaces` in the repository field.
5. Select the `mapping` template.
6. Select `2 cores, 8GB RAM, 32 GB storage` option.

It should take a few minutes to create the Codespace. Once it's ready, you will be able to access the terminal and files directly from VS Code.

If you need a graphical interface, you can use TigerVNC to connect to the Codespace. Open TigerVNC and enter the following address: `localhost:5901`. Once you've entered the password, you should see the full desktop environment of your Codespace.