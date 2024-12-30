#!/bin/bash
# This script reads D, Dstar and F files from each algorithm and transform them to T1w space.

# Define input arguments
csv_file=$1
reference_image=$2
xfm_file=$3
interpolation_method=$4
output_path=$5
output_csv=$6

echo "Algorithm,Metric" >> "$output_csv"

# Read the CSV file column by column
column_count=$(head -n 1 "$csv_file" | tr ',' '\n' | wc -l)
for ((i = 1; i <= column_count; i++))
do
    column_name=$(awk -F',' -v col="$i" 'NR==1 {print $col}' "$csv_file")
    image_paths=$(awk -F',' -v col="$i" 'NR>1 {print $col}' "$csv_file")

    # Process each image path in the column
    while read -r image_path
    do

        # Stripping just the metric name.
        filename=$(echo "$image_path" | awk -F'concatted_' '{print $2}')

        echo "Registering ${filename} from ${column_name} algorithm to T1w space"

        transformed_filename="${output_path}_${column_name}/transformed_${filename}"
        # Run antsApplyTransforms for the current image
        antsApplyTransforms -d 3 --input-image-type 3 -i "$image_path" -r "$reference_image" -t "$xfm_file" -o "$transformed_filename" -n "$interpolation_method"

        # Save the column name and transformed filename to the output CSV file
        echo "$column_name,$transformed_filename" >> "$output_csv"

    done <<< "$image_paths"
done


