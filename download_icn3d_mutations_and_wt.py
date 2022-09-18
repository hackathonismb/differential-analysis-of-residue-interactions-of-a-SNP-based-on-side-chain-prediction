"""Download PDB and Mutation files from iCn3D"""

import argparse
import csv
import os
import re
import time
import sys

from selenium import webdriver
from selenium.webdriver.firefox.options import Options


CRAWL_DELAY = 5
MAX_SLEEP = 20
MAX_ATTEMPTS = 2

def argument_parser():
    """Parse arguments for the distance_to_features scipt."""
    parser = argparse.ArgumentParser()
    required_arguments = parser.add_argument_group("Required Arguments")
    required_arguments.add_argument("-i",
                                    "--input_file",
                                    required=True,
                                    type=str,
                                    help="Input file with PDB mutation information.")
    required_arguments.add_argument("-o",
                                    "--output_folder",
                                    required=True,
                                    type=str,
                                    help="Output folder for downloads.")
    parser.add_argument("-v",
                        "--verbose",
                        action="store_true",
                        help="Prints updates on the download progress.")
    args = parser.parse_args()
    return args


def read_csv_file(filename, skip_header=True):
    """Reads a csv file as a nested list."""
    with open(filename, newline='') as csvfile:
        reader = csv.reader(csvfile)
        if skip_header:
            next(reader) # skip header
        rows = list(reader)
        return rows


def download_file_from_icn3d(api_call,
                             expected_file_path,
                             options,
                             crawl_delay=CRAWL_DELAY,
                             max_sleep=MAX_SLEEP,
                             max_attempts=MAX_ATTEMPTS):
    """Downloads a file from icn3d and checks download completion.

    Given an API call to iCn3D as a web address string, the expected
    file location, options for the webdriver, and the crawl delay,
    downloads a file from iCn3D and checks if the file was downloaded.
    If not, returns a FileNotFoundError.
    """
    attempt = 0
    is_downloaded = False
    # Loop until file is downloaded or max attempts reached
    while attempt < max_attempts and is_downloaded is False:
        # Open webdriver
        with webdriver.Firefox(options=options) as driver:
            # Open icn3d with api call and wait
            driver.get(api_call)
            time.sleep(crawl_delay)
            # Update download and sleep status
            is_downloaded = os.path.isfile(expected_file_path)
            total_sleep = crawl_delay
            # If the file has not downloaded yet, wait
            while is_downloaded is False and total_sleep <= max_sleep:
                time.sleep(crawl_delay)
                total_sleep += crawl_delay
                is_downloaded = os.path.isfile(expected_file_path)
        attempt += 1
    if is_downloaded is False:
        expected_file = expected_file_path.split(r'/')[-1]
        raise FileNotFoundError(f"{expected_file} failed to download after {attempt} attempt(s).")

def set_webdriver_options(download_folder, headless=True):
    """Set the driver options for the firefox selenium driver.

    Sets the headless and download location options for the
    firefox selenium driver. By default, it will run headless and
    use the default download folder.
    """
    options = Options()
    # Sets the download folder to a specified path
    options.set_preference("browser.download.folderList", 2)
    options.set_preference("browser.download.dir", download_folder)
    # If headless, adds a --headless argument to the Options object
    if headless:
        options.add_argument( '--headless' )
    return options


def print_progress(expected_file, idx, total):
    print(f"Downloading {expected_file}. File {idx + 1}/{total}   ",
          end='\r')


def download_mutations(rows,
                       webdriver_options,
                       verbose=False):
    for idx, row in enumerate(rows):
        pdb_id = row[1]
        chain = row[2]
        residue = re.search(r"[0-9]+", row[3]).group(0)
        mutation = row[3][-1]
        mutation_api = ("https://www.ncbi.nlm.nih.gov/Structure/icn3d/full.html?"
                        f"pdbid={pdb_id}&command=scap%20pdb%20{pdb_id}_{chain}_{residue}_{mutation};"
                        "%20select%20displayed%20set")
        # File downloaded, not correct format?
        #mutation_api = ("https://www.ncbi.nlm.nih.gov/Structure/icn3d/full.html?"
        #               f"pdbid={pdb_id}&command=scap%20pdb%20{pdb_id}_{chain}_{residue}_{mutation}")
        expected_file = f"{pdb_id}_{chain}_{residue}_{mutation}.pdb"
        download_folder = webdriver_options._preferences["browser.download.dir"]
        expected_file_path = f"{download_folder}/{expected_file}"
        if verbose:
            print_progress(expected_file, idx, len(rows))
        # If the file has not been downloaded before, download it
        if os.path.isfile(expected_file_path) is False:
            download_file_from_icn3d(mutation_api,
                                     expected_file_path,
                                     webdriver_options)
    if verbose:
        print("\nDownload of PDB mutation file(s) complete.")


def download_pdb_files(pdb_ids,
                       webdriver_options,
                       verbose=False):
    for idx, pdb_id in enumerate(pdb_ids):
        wt_api = ("https://www.ncbi.nlm.nih.gov/Structure/icn3d/full.html?"
                  f"pdbid={pdb_id}&command=export%20pdb")
        expected_file = f"{pdb_id}_icn3d.pdb"
        download_folder = webdriver_options._preferences["browser.download.dir"]
        expected_file_path = f"{download_folder}/{expected_file}"
        if verbose:
            print_progress(expected_file, idx, len(pdb_ids))
        if os.path.isfile(expected_file_path) is False:
            download_file_from_icn3d(wt_api,
                                     expected_file_path,
                                     webdriver_options)
    if verbose:
        print("\nDownload of PDB file(s) complete.")


def main():
    args = argument_parser()
    output_folder_abs = os.path.abspath(args.output_folder)
    options = set_webdriver_options(output_folder_abs)
    mutation_rows = read_csv_file(args.input_file)
    download_mutations(mutation_rows, options, verbose=args.verbose)
    unique_pdb_ids = list(set(map(lambda x: x[1], mutation_rows)))
    download_pdb_files(unique_pdb_ids, options, verbose=args.verbose)


if __name__ == "__main__":
    main()

