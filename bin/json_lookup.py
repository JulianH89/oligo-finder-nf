#!/usr/bin/env python

import json
import argparse
import sys

def get_accessions(json_file_path, data_id, mismatch_level):
    """
    Looks up accessions in a JSON file based on an ID and mismatch level.

    Args:
        json_file_path (str): The path to the JSON file.
        data_id (str): The ID to search for (e.g., '16').
        mismatch_level (str): The mismatch level to search for (e.g., '0').

    Returns:
        list: A list of accession strings if found, otherwise None.
    """
    try:
        with open(json_file_path, 'r') as f:
            data = json.load(f)
        
        # Navigate through the nested dictionary to find the accessions
        return data[data_id]['mismatch_level'][mismatch_level]['accessions']

    except FileNotFoundError:
        print(f"Error: The file '{json_file_path}' was not found.", file=sys.stderr)
        return None
    except KeyError:
        print(f"Error: Could not find the path for ID '{data_id}' and mismatch level '{mismatch_level}'.", file=sys.stderr)
        return None
    except json.JSONDecodeError:
        print(f"Error: The file '{json_file_path}' is not a valid JSON file.", file=sys.stderr)
        return None

if __name__ == "__main__":
    # Set up the command-line argument parser
    parser = argparse.ArgumentParser(
        description="A simple script to look up accessions in a JSON file."
    )
    parser.add_argument("--json", help="Path to the input JSON file.")
    parser.add_argument("--id", help="The ID to look up (e.g., 16).")
    parser.add_argument("--mismatch_level", help="The mismatch level (e.g., 0).")

    args = parser.parse_args()

    # Call the function with command-line arguments
    result = get_accessions(args.json, args.id, args.mismatch_level)

    # Print the results to the console
    if result:
        for accession in result:
            print(accession)
    else:
        # Error messages are printed from within the function
        sys.exit(1)