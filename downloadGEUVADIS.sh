

mkdir rawReads
cd rawReads
# Parse and download files from manifest
while IFS=$'\t' read -r url md5 project type molecule sample population; do
    filename=$(basename "$url")
    
    echo "Downloading: $filename"
    wget "$url" -O "$filename"
    
    # Verify MD5 checksum
    echo "Verifying checksum..."
    calculated_md5=$(md5sum "$filename" | cut -d' ' -f1)
    
    if [ "$calculated_md5" = "$md5" ]; then
        echo "✓ Checksum verified for $filename"
    else
        echo "✗ Checksum mismatch for $filename"
        echo "Expected: $md5"
        echo "Got: $calculated_md5"
    fi
done < /scratch/vasccell/cs806/geuvadis/geuvadis_Scripts/igsr_Geuvadis.tsv.tsv

