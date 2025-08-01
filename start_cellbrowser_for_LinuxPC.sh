#!/bin/bash
# filepath: /vol/research/scratch1/NOBACKUP/hh01116/10X_DIPG/start_cellbrowser.sh

set -e  # Exit on any error

echo "=== Cell Browser Setup and Launch Script ==="

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to check if cellbrowser is installed
check_cellbrowser() {
    if command_exists cbUpgrade; then
        echo "✓ cellbrowser is already installed"
        return 0
    else
        echo "✗ cellbrowser is not installed"
        return 1
    fi
}

# Function to ask for user confirmation
ask_confirmation() {
    local message="$1"
    echo "$message"
    read -p "Do you want to proceed? (y/N): " -n 1 -r
    echo    # Move to a new line
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        return 0
    else
        echo "Installation cancelled by user."
        return 1
    fi
}

# Function to install cellbrowser via conda
install_via_conda() {
    echo "Attempting to install cellbrowser via conda..."
    if ! ask_confirmation "This will install ucsc-cell-browser from bioconda channel."; then
        return 1
    fi
    conda install -c bioconda ucsc-cell-browser -y
    if [ $? -eq 0 ]; then
        echo "✓ Successfully installed cellbrowser via conda"
        return 0
    else
        echo "✗ Failed to install cellbrowser via conda"
        return 1
    fi
}

# Function to install cellbrowser via pip
install_via_pip() {
    echo "Attempting to install cellbrowser via pip..."
    if ! ask_confirmation "This will install cellbrowser via pip (--user flag will be used)."; then
        return 1
    fi
    pip install --user cellbrowser
    if [ $? -eq 0 ]; then
        echo "✓ Successfully installed cellbrowser via pip"
        # Add user's local bin to PATH if not already there
        if [[ ":$PATH:" != *":$HOME/.local/bin:"* ]]; then
            export PATH="$HOME/.local/bin:$PATH"
            echo "Added $HOME/.local/bin to PATH"
        fi
        return 0
    else
        echo "✗ Failed to install cellbrowser via pip"
        return 1
    fi
}

# Main installation logic
install_cellbrowser() {
    echo "Installing cellbrowser..."
    
    # Try conda first if available
    if command_exists conda; then
        echo "Found conda, trying conda installation first..."
        if install_via_conda; then
            return 0
        fi
    fi
    
    # Try pip if conda failed or not available
    if command_exists pip; then
        echo "Trying pip installation..."
        if install_via_pip; then
            return 0
        fi
    fi
    
    echo "✗ Failed to install cellbrowser via both conda and pip"
    return 1
}

# Check for package managers
echo "Checking for package managers..."
has_conda=false
has_pip=false

if command_exists conda; then
    echo "✓ conda is available"
    has_conda=true
else
    echo "✗ conda is not available"
fi

if command_exists pip; then
    echo "✓ pip is available"
    has_pip=true
else
    echo "✗ pip is not available"
fi

# Exit if neither package manager is available
if [ "$has_conda" = false ] && [ "$has_pip" = false ]; then
    echo "Error: Neither conda nor pip is available. Please install one of them first."
    exit 1
fi

# Check if cellbrowser is installed, install if not
if ! check_cellbrowser; then
    echo "cellbrowser not found, attempting installation..."
    if ! install_cellbrowser; then
        echo "Error: Failed to install cellbrowser"
        exit 1
    fi
fi

# Verify installation
if ! check_cellbrowser; then
    echo "Error: cellbrowser installation verification failed"
    exit 1
fi

# Check if cellbrowser directory exists
CELLBROWSER_DIR="./output/html_reports/cellbrowser"
if [ ! -d "$CELLBROWSER_DIR" ]; then
    echo "Warning: cellbrowser directory '$CELLBROWSER_DIR' does not exist"
    echo "Please make sure you have exported data using the Cell Browser exporter first"
    echo "Looking for alternative directories..."
    
    # Look for alternative cellbrowser directories
    ALT_DIRS=("./html_reports/cellbrowser" "./cellbrowser")
    
    for dir in "${ALT_DIRS[@]}"; do
        if [ -d "$dir" ]; then
            echo "Found alternative directory: $dir"
            CELLBROWSER_DIR="$dir"
            break
        fi
    done
    
    if [ ! -d "$CELLBROWSER_DIR" ]; then
        echo "Error: No cellbrowser directory found. Please export data first."
        exit 1
    fi
fi

echo "Using cellbrowser directory: $CELLBROWSER_DIR"

# Start the cellbrowser server
echo "Starting Cell Browser server..."
echo "Server will be available at: http://localhost:8899"
echo "Press Ctrl+C to stop the server"
echo ""

cbUpgrade -o "$CELLBROWSER_DIR" -p 8899