# symbol_response = perform_request("lookup/symbol/homo_sapiens/Q8WZ42")
# id_ = "ENSG00000155657"
# id_ = "ENST00000342992"
# lookup_response = perform_request(f"lookup/id/{id_}")
#
# cdna_sequence_response = perform_request(f"sequence/id/{id_}?type=cdna")
#
# cds_sequence_response = perform_request(f"sequence/id/{id_}?type=cds")

# f = open("ttn-202.txt", "r")
# content = f.read()
# raw = content.rstrip().replace("\n", "")
#
# f = open("uniprot.txt", "r")
# uniprot_content = f.read().rstrip().replace("\n", "")
#
# fr = open("ttn-202-result.txt", "r")
# content2 = fr.read()

# region = f"{lookup_response['start']}..{lookup_response['end']}"
# cdna_region = perform_request(f"map/cdna/{id_}/{region}")
#
# translation = perform_request(f"map/translation/{id_}/{region}") # Tira 500
#
# translation