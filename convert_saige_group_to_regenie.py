import argparse
import gzip

def convert_annotation_string(input_str):
    # Split the string by ':' and iterate over each part
    annotations = input_str.split(',')
    output_lines = []

    for annotation in annotations:

        parts = annotation.split(':')
        output_lines.append(f"{annotation} {','.join(parts)}")

    return output_lines

def process_file(input_file, annotation_file, set_list_file, annotation_string):
    gene_variants = {}
    annotations = []

    # Convert the colon-separated string and write to a file
    converted_lines = convert_annotation_string(annotation_string)
    with open('mask.txt', 'w') as conv_file:
        for line in converted_lines:
            conv_file.write(line + '\n')

    with gzip.open(input_file, 'r') as file:
        for line in file:
            parts = line.strip().split()
            gene_name = parts[0]
            var_info = parts[1] == 'var'
            anno_info = parts[1] == 'anno'

            if var_info:
                # Process variant information
                variants = parts[2:]
                gene_variants[gene_name] = variants
            elif anno_info:
                # Process annotation information
                for anno, variant in zip(parts[2:], gene_variants.get(gene_name, [])):
                    chrom, pos, ref, alt = variant.split(':')
                    annotations.append(f"{chrom}:{pos}:{ref}:{alt} {gene_name} {anno}")

    # Writing to Annotation File
    with open(annotation_file, 'w') as anno_file:
        for annotation in annotations:
            anno_file.write(annotation + '\n')

    # Writing to Set List File
    with open(set_list_file, 'w') as set_file:
        for gene, variants in gene_variants.items():
            set_line = f"{gene} {' '.join(variants[0].split(':')[0:2])} {','.join(variants)}"
            set_file.write(set_line + '\n')

def main():
    parser = argparse.ArgumentParser(description='Convert gene data file to annotation and set list files.')
    parser.add_argument('input_file', help='Input file path')
    parser.add_argument('annotation_file', help='Output file path for annotations')
    parser.add_argument('set_list_file', help='Output file path for set list')
    parser.add_argument('annotation_string', help='String-separated string to be converted')

    args = parser.parse_args()
    process_file(args.input_file, args.annotation_file, args.set_list_file, args.annotation_string)

if __name__ == "__main__":
    main()
