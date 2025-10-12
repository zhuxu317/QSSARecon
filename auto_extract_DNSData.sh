#!/bin/bash
target="/data/ZhuXu/Cantera/NeuralRecon/cases/H2_DNS/DNSData.tar.gz"
prev_size=0

echo "🔍 Monitoring file growth..."
while true; do
  size=$(stat -c%s "$target" 2>/dev/null || echo 0)
  echo "$(date '+%H:%M:%S')  size = $size bytes"
  if [ "$size" -eq "$prev_size" ] && [ "$size" -ne 0 ]; then
    echo "✅ Size stable — starting extraction..."
    cd /data/ZhuXu/Cantera/NeuralRecon/cases/H2_DNS || exit 1
    if command -v pigz >/dev/null 2>&1; then
      tar -I pigz -xvf DNSData.tar.gz
    else
      tar -xvzf DNSData.tar.gz
    fi
    echo "✅ Extraction done."
    break
  fi
  prev_size=$size
  sleep 60
done