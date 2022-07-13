from uniprot_ID_mapping import *
import requests
import json
import argparse

def argument_parser():
    parser = argparse.ArgumentParser()
    required_arguments = parser.add_argument_group("Required Arguments")
    required_arguments.add_argument("id",
                                    type=str,
                                    help="query ID")
    required_arguments.add_argument("db",
                                    type=str,
                                    help="data base can be either uniprot or pdb")
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        help="output json filename")
    args = parser.parse_args()
    return args


def query_UniProt(uniprotkb):
    # Query UniProt through the REST interface
    API_URL = f'https://rest.uniprot.org/uniprot/{uniprotkb}'
    return requests.get(API_URL).json()


def retrieve_feature_regions(data):
    features_dict = {}
    for feature in data['features']:
        if feature['type'] in ['Binding site',
                               'Region',
                               'Disulfide bond',
                               'Glycosylation',
                               'Lipidation']:
            if feature['type'] not in features_dict:
                features_dict[feature['type']] = []
            start = feature['location']['start']['value']
            end = feature['location']['end']['value']
            if start != end:
                if feature['type'] == 'Disulfide bond':
                    features_dict[feature['type']].append((start, end))
                else:
                    features_dict[feature['type']].append([start, end])
            else:
                # those annotations with a single residue are set to have
                # the same "start" and "end"
                features_dict[feature['type']].append(start)
    return features_dict


def main():
    args = argument_parser()
    if args.db == 'uniprot':
        data = query_UniProt(args.id)
        if 'messages' in data:
            raise ValueError(f"invalid uniprotKB: {args.id}")
        features_dict = retrieve_feature_regions(data)
    elif args.db == 'pdb':
        job_id = submit_id_mapping(from_db="PDB", to_db="UniProtKB", ids=[args.id])
        if check_id_mapping_results_ready(job_id):
            link = get_id_mapping_results_link(job_id)
            results = get_id_mapping_results_search(link)
        features_dict = {}
        for hit in results['results']:
            if hit['from'] not in features_dict:
                features_dict[hit['from']] = {}
            features_dict[hit['from']][hit['to']['primaryAccession']] =  retrieve_feature_regions(hit['to'])
    else:
        raise ValueError(f"invalid db {args.db}. Choose between uniprot and pdb")
    if args.output != None and '.json' in args.output:
        filename = args.output
    else: 
        filename = args.id + '_regions.json'
    with open(filename, "w") as f:
        f.write(json.dumps(features_dict))


if __name__ == "__main__":
    main()