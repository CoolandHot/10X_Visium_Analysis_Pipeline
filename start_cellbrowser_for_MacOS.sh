#!/bin/bash

# ============================================================================
#  Cellbrowser Upgrade Script for macOS
# ============================================================================
#  This script automates running the cbUpgrade command from the cellbrowser 
#  package. It tries multiple methods to find cbUpgrade and executes it.
#  Includes macOS-specific paths like Homebrew and user Library locations.
# ============================================================================

# ============================================================================
#  HELPER FUNCTIONS
# ============================================================================

find_via_python_location() {
    # Attempt to find cbUpgrade via Python executable location
    local python_exe
    local python_dir
    local bin_dir
    
    if python_exe=$(python3 -c "import sys; print(sys.executable)" 2>/dev/null) || 
       python_exe=$(python -c "import sys; print(sys.executable)" 2>/dev/null); then
        python_dir=$(dirname "$python_exe")
        
        # Try bin directory (common in virtual environments and user installs)
        bin_dir="$python_dir/bin"
        if [[ -f "$bin_dir/cbUpgrade" ]]; then
            echo "      SUCCESS: Found at $bin_dir/cbUpgrade"
            CB_UPGRADE_CMD="$bin_dir/cbUpgrade"
            return
        fi
        
        # Try ../bin directory (system installations)
        bin_dir="$(dirname "$python_dir")/bin"
        if [[ -f "$bin_dir/cbUpgrade" ]]; then
            echo "      SUCCESS: Found at $bin_dir/cbUpgrade"
            CB_UPGRADE_CMD="$bin_dir/cbUpgrade"
            return
        fi
    fi
}

find_via_macos_locations() {
    # Check macOS-specific installation locations
    local locations=(
        "/opt/homebrew/bin/cbUpgrade"           # Apple Silicon Homebrew
        "/usr/local/bin/cbUpgrade"              # Intel Homebrew
        "$HOME/Library/Python/*/bin/cbUpgrade" # User Python packages
        "$HOME/.local/bin/cbUpgrade"            # User local installs
        "/Library/Frameworks/Python.framework/Versions/*/bin/cbUpgrade" # Python.org installs
    )
    
    for location in "${locations[@]}"; do
        # Handle glob patterns
        if [[ "$location" == *"*"* ]]; then
            for expanded in $location; do
                if [[ -f "$expanded" ]]; then
                    echo "      SUCCESS: Found at $expanded"
                    CB_UPGRADE_CMD="$expanded"
                    return
                fi
            done
        else
            if [[ -f "$location" ]]; then
                echo "      SUCCESS: Found at $location"
                CB_UPGRADE_CMD="$location"
                return
            fi
        fi
    done
}

find_via_pip_show() {
    # Attempt to find cbUpgrade via pip package information
    local temp_file="/tmp/pip_show_$$.txt"
    local pkg_location
    local python_path
    local bin_dir
    
    # Try both pip3 and pip
    if pip3 show cellbrowser > "$temp_file" 2>&1 || pip show cellbrowser > "$temp_file" 2>&1; then
        # Check if output file has content
        if [[ -s "$temp_file" ]]; then
            # Extract package location
            pkg_location=$(grep "Location:" "$temp_file" | cut -d' ' -f2-)
            
            if [[ -n "$pkg_location" ]]; then
                # Derive Python bin directory from package location
                # Handle common macOS patterns
                if [[ "$pkg_location" == */site-packages ]]; then
                    python_path="${pkg_location%/site-packages}"
                elif [[ "$pkg_location" == */lib/python*/site-packages ]]; then
                    python_path="${pkg_location%/lib/python*/site-packages}"
                elif [[ "$pkg_location" == *Library/Python* ]]; then
                    # Handle ~/Library/Python/3.x/lib/python/site-packages
                    python_path=$(echo "$pkg_location" | sed 's|/lib/python/site-packages||')
                else
                    python_path="$pkg_location"
                fi
                
                # Try bin directory
                bin_dir="$python_path/bin"
                if [[ -f "$bin_dir/cbUpgrade" ]]; then
                    echo "      SUCCESS: Found at $bin_dir/cbUpgrade"
                    CB_UPGRADE_CMD="$bin_dir/cbUpgrade"
                fi
            fi
        fi
    fi
    
    [[ -f "$temp_file" ]] && rm -f "$temp_file"
}

exec_upgrade() {
    echo
    echo "============================================================================"
    echo "  Cellbrowser Tool Located: $CB_UPGRADE_CMD"
    echo "============================================================================"

    # Determine script directory
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

    # Try to auto-detect input directory
    AUTO_INPUT=""
    if [[ -d "$SCRIPT_DIR/cellbrowser" ]]; then
        read -p "Detected 'cellbrowser' subfolder in script directory. Use '$SCRIPT_DIR/cellbrowser' as input directory? [Y/n]: " yn
        case "$yn" in
            [Nn]*) ;;
            *) AUTO_INPUT="$SCRIPT_DIR/cellbrowser";;
        esac
    elif [[ -f "$SCRIPT_DIR/index.html" ]]; then
        read -p "Detected 'index.html' in script directory. Use '$SCRIPT_DIR' as input directory? [Y/n]: " yn
        case "$yn" in
            [Nn]*) ;;
            *) AUTO_INPUT="$SCRIPT_DIR";;
        esac
    fi

    if [[ -n "$AUTO_INPUT" ]]; then
        INPUT_DIR="$AUTO_INPUT"
    else
        # Get input directory from user
        read -p "Enter the cellbrowser input directory: " INPUT_DIR
    fi

    # Validate directory exists
    if [[ ! -d "$INPUT_DIR" ]]; then
        echo
        echo "ERROR: Directory '$INPUT_DIR' does not exist."
        echo "Please check the path and try again."
        cleanup_and_exit
    fi

    # Execute the upgrade command
    PORT="8899"
    echo
    echo "------------------------------------------------------------"
    echo "Executing: $CB_UPGRADE_CMD -o \"$INPUT_DIR\" -p $PORT"
    echo "------------------------------------------------------------"
    echo

    $CB_UPGRADE_CMD -o "$INPUT_DIR" -p $PORT
    EXIT_CODE=$?

    # Check execution result
    if [[ $EXIT_CODE -ne 0 ]]; then
        echo
        echo "WARNING: Command completed with error code $EXIT_CODE"
    else
        echo
        echo "SUCCESS: Command completed successfully."
        echo "The cellbrowser should be accessible at: http://localhost:$PORT"
    fi
    cleanup_and_exit
}

cleanup_and_exit() {
    # Clean up any temporary files
    rm -f /tmp/python_path_$$.txt /tmp/pip_show_$$.txt
    
    echo
    echo "============================================================================"
    echo "  Script finished."
    echo "============================================================================"
    read -p "Press Enter to continue..."
    exit 0
}

# ============================================================================
#  MAIN SCRIPT EXECUTION
# ============================================================================

echo "============================================================================"
echo "  Starting Cellbrowser Upgrade Script for macOS"
echo "============================================================================"
echo

CB_UPGRADE_CMD=""

# ----------------------------------------------------------------------------
#  METHOD 1: Check if cbUpgrade is directly available in PATH
# ----------------------------------------------------------------------------
echo "[1/5] Checking if cbUpgrade is available in PATH..."
if command -v cbUpgrade >/dev/null 2>&1; then
    echo "      SUCCESS: Found cbUpgrade in PATH"
    CB_UPGRADE_CMD="cbUpgrade"
    exec_upgrade
else
    echo "      Not found in PATH"
fi

# ----------------------------------------------------------------------------
#  METHOD 2: Check macOS-specific installation locations
# ----------------------------------------------------------------------------
echo "[2/5] Checking macOS-specific locations (Homebrew, user Library)..."
find_via_macos_locations
if [[ -n "$CB_UPGRADE_CMD" ]]; then
    exec_upgrade
else
    echo "      Not found in macOS-specific locations"
fi

# ----------------------------------------------------------------------------
#  METHOD 3: Locate cbUpgrade via Python installation directory
# ----------------------------------------------------------------------------
echo "[3/5] Checking Python installation directory..."
find_via_python_location
if [[ -n "$CB_UPGRADE_CMD" ]]; then
    exec_upgrade
else
    echo "      Not found via Python location"
fi

# ----------------------------------------------------------------------------
#  METHOD 4: Locate cbUpgrade via pip package information
# ----------------------------------------------------------------------------
echo "[4/5] Checking cellbrowser package location..."
find_via_pip_show
if [[ -n "$CB_UPGRADE_CMD" ]]; then
    exec_upgrade
else
    echo "      Not found via pip show"
fi

# ----------------------------------------------------------------------------
#  METHOD 5: Try running cellbrowser as Python module
# ----------------------------------------------------------------------------
echo "[5/5] Trying cellbrowser as Python module..."
if python3 -m cellbrowser.cbUpgrade --help >/dev/null 2>&1; then
    echo "      SUCCESS: Found cellbrowser module (python3)"
    CB_UPGRADE_CMD="python3 -m cellbrowser.cbUpgrade"
    exec_upgrade
elif python -m cellbrowser.cbUpgrade --help >/dev/null 2>&1; then
    echo "      SUCCESS: Found cellbrowser module (python)"
    CB_UPGRADE_CMD="python -m cellbrowser.cbUpgrade"
    exec_upgrade
else
    echo "      Module not available"
fi

# ----------------------------------------------------------------------------
#  ERROR: All methods failed
# ----------------------------------------------------------------------------
echo
echo "============================================================================"
echo "  ERROR: Could not locate cbUpgrade executable"
echo "============================================================================"
echo "  Please ensure:"
echo "    1. Python is installed (try: brew install python)"
echo "    2. cellbrowser package is installed: pip3 install cellbrowser"
echo "    3. Check if Homebrew is installed: /bin/bash -c \"\$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)\""
echo "    4. Python bin directory is in your PATH"
echo "============================================================================"
cleanup_and_exit
