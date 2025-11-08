#!/usr/bin/env bash
#SBATCH --job-name=install_relate
#SBATCH --output=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/install_relate_%j.out
#SBATCH --error=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/install_relate_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

# Install Relate tool for selection analysis verification

SOFTWARE_DIR="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/software"
RELATE_DIR="$SOFTWARE_DIR/relate"

mkdir -p "$SOFTWARE_DIR"
cd "$SOFTWARE_DIR"

echo "=========================================="
echo "Installing Relate"
echo "=========================================="

# Download Relate (latest version)
RELATE_VERSION="v1.2.1"
RELATE_URL="https://github.com/MyersGroup/relate/releases/download/${RELATE_VERSION}/relate_${RELATE_VERSION}_x86_64_static.tgz"

echo "Downloading Relate ${RELATE_VERSION}..."
wget "$RELATE_URL" -O relate.tgz

if [[ $? -ne 0 ]]; then
    echo "ERROR: Failed to download Relate"
    exit 1
fi

echo "Extracting Relate..."
tar -xzf relate.tgz
rm relate.tgz

# Rename directory
mv relate_*_x86_64_static "$RELATE_DIR" 2>/dev/null || echo "Directory already exists"

echo ""
echo "✓ Relate installed successfully"
echo ""
echo "Installation directory: $RELATE_DIR"
echo ""
echo "Add to your PATH:"
echo "  export PATH=\"$RELATE_DIR/bin:\$PATH\""
echo ""
echo "Or add to your ~/.bashrc:"
echo "  echo 'export PATH=\"$RELATE_DIR/bin:\$PATH\"' >> ~/.bashrc"
echo ""
echo "Test installation:"
echo "  $RELATE_DIR/bin/RelateFileFormats --help"
