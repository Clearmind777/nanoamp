#!/bin/sh
# Install a small `nanoamp` wrapper into a directory on PATH.
set -e

DEST="${1:-$HOME/.local/bin}"
mkdir -p "$DEST"

cat > "$DEST/nanoamp" <<'EOF'
#!/bin/sh
exec Rscript --vanilla -e "library(nanoamp); nanoamp_cli()" "$@"
EOF
chmod +x "$DEST/nanoamp"

echo "Installed: $DEST/nanoamp"
case ":$PATH:" in
  *":$DEST:"*) ;;
  *) echo "Note: add $DEST to PATH, for example: export PATH=\"$DEST:\$PATH\"" ;;
esac
