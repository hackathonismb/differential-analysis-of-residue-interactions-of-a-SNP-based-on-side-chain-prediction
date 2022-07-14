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
                               'Lipidation',
                               'Metal binding',
                               'Active site',
                               'Cross-link']:
            if feature['type'] not in features_dict:
                features_dict[feature['type']] = []
            start = feature['location']['start']['value']
            end = feature['location']['end']['value']
            if start != end:
                if feature['type'] == 'Disulfide bond':
                    coord = (start, end)
                else:
                    coord = [start, end]
                
                if coord not in features_dict[feature['type']]:
                    features_dict[feature['type']].append(coord)
            else:
                # those annotations with a single residue are set to have
                # the same "start" and "end"
                if start not in features_dict[feature['type']]:
                    features_dict[feature['type']].append(start)
    return features_dict


def main():
    args = argument_parser()
    with open('Combined_PTMs.json', "r") as f:
        combined_PTMs = json.load(f) 
    if args.db == 'uniprot':
        data = query_UniProt(args.id)
        if 'messages' in data:
            raise ValueError(f"invalid uniprotKB: {args.id}")
        features_dict = retrieve_feature_regions(data)
        if args.id in combined_PTMs:
            features_dict = {**features_dict, **combined_PTMs[args.id]}
    elif args.db == 'pdb':
        job_id = submit_id_mapping(from_db="PDB", to_db="UniProtKB", ids=[args.id])
        if check_id_mapping_results_ready(job_id):
            link = get_id_mapping_results_link(job_id)
            results = get_id_mapping_results_search(link)
        features_dict = {}
        for hit in results['results']:
            id_ = hit['to']['primaryAccession']
            feature_hits_dict = retrieve_feature_regions(hit['to'])
            if id_ in combined_PTMs:
                feature_hits_dict = {**feature_hits_dict, **combined_PTMs[id_]}
            features_dict[hit['to']['primaryAccession']] =  feature_hits_dict 
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