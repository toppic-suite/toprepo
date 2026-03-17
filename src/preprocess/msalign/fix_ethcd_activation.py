#!/usr/bin/env python3
"""
Script to replace "ACTIVATION=ETD" with "ACTIVATION=ETHCD" in files
listed in a specified file list.

Usage: python fix_ethcd_activation.py <file_list>
"""

import os
import sys

def process_file(filepath):
    """Replace ACTIVATION=ETD with ACTIVATION=ETHCD in a file."""
    try:
        # Read the file content
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()

        # Count occurrences before replacement
        count = content.count('ACTIVATION=ETD')

        if count == 0:
            print(f"  No matches found in {filepath}")
            return 0

        # Replace the text
        new_content = content.replace('ACTIVATION=ETD', 'ACTIVATION=ETHCD')

        # Write back to the file
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(new_content)

        print(f"  Replaced {count} occurrence(s) in {filepath}")
        return count

    except FileNotFoundError:
        print(f"  ERROR: File not found: {filepath}")
        return 0
    except Exception as e:
        print(f"  ERROR processing {filepath}: {e}")
        return 0

def main():
    if len(sys.argv) != 2:
        print("Usage: python fix_ethcd_activation.py <file_list>")
        sys.exit(1)

    file_list_path = sys.argv[1]

    if not os.path.exists(file_list_path):
        print(f"ERROR: {file_list_path} not found!")
        sys.exit(1)

    # Read the list of files
    with open(file_list_path, 'r') as f:
        files = [line.strip() for line in f if line.strip()]

    print(f"Processing {len(files)} files...\n")

    total_replacements = 0
    processed_count = 0

    for filepath in files:
        replacements = process_file(filepath)
        total_replacements += replacements
        if replacements > 0:
            processed_count += 1

    print(f"\nSummary:")
    print(f"  Files processed: {processed_count}/{len(files)}")
    print(f"  Total replacements: {total_replacements}")

if __name__ == '__main__':
    main()
