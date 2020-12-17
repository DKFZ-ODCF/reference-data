#!/usr/bin/env python3

import requests, sys, json


def request_from_api(ext, query):
    server = "https://rest.ensembl.org"

    response = requests.get(server + ext + query, headers={"Content-Type": "application/json"})

    if not response.ok:
        response.raise_for_status()
        sys.exit()

    return response


def arg_to_str(args_list):
    args_as_str = " ".join(args_list)
    args_as_str = args_as_str.replace(" ", "%20")
    return args_as_str


def display_taxonomy_result(taxdict):
    print("Select your assembly by the number:")
    for i, elem in enumerate(taxdict):
        print("{i} - {dispname}: {assemblyname}".format(i=i,
                                                        dispname=elem.get("display_name"),
                                                        assemblyname=elem.get("name")))
    return True


def request_user_selection(taxdict):
    while True:
        try:
            select = int(input("Enter number: "))
            if select not in range(len(taxdict)):
                print("Number not in range")
                continue
        except ValueError:
            print("Sorry, I didn't understand that.")
            # better try again... Return to the start of the loop
            continue
        else:
            break

    assembly_name = taxdict[select].get("name")
    print(f"You selected {select}: {assembly_name}")
    assembly_json_path = f"resources/{assembly_name}.json"

    r = request_from_api("/info/genomes/", assembly_name)
    with open(assembly_json_path, "w") as json_out:
        json.dump(r.json(), fp=json_out, indent=4)

    print(f"The JSON file for the selected assembly has been stored to {assembly_json_path} ."
          "If you chose to run the pipeline with it, don't forget to commit the file to the repository,"
          " and copy it 'resouces/current.json'")


species = arg_to_str(sys.argv[1:])
r = request_from_api("/info/genomes/taxonomy/", species)
species_response_dict = r.json()

display_taxonomy_result(species_response_dict)
request_user_selection(species_response_dict)
