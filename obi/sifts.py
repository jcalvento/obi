import gzip
import os
import shutil
from contextlib import closing
from urllib import request
from xml.etree import ElementTree

from obi.utils import root_path


def _create_results_dir(pdb_root_id):
    dir_path = f"{root_path()}/sifts"
    if not os.path.isdir(dir_path):
        os.mkdir(dir_path)
    pdb_dir_path = f"{dir_path}/{pdb_root_id}"
    if not os.path.isdir(pdb_dir_path):
        os.mkdir(pdb_dir_path)
    return dir_path


class Sifts:
    def map_to(self, pdb_id):
        pdb_xml = self.__download_xml(pdb_id)
        self.__parse_sifts(pdb_xml)

    def __download_xml(self, pdb_id):
        pdb_root_id = pdb_id[1:3]
        pdb_xml_path = f"{pdb_root_id}/{pdb_id}.xml"
        xml_file = f"{root_path()}/sifts/{pdb_xml_path}"
        if not os.path.isfile(xml_file):
            print("Downloading pdb xml")
            _create_results_dir(pdb_root_id)
            xml_url = f'http://ftp.ebi.ac.uk/pub/databases/msd/sifts/split_xml/{pdb_xml_path}.gz'
            print(f"Fetching xml from {xml_url}")
            try:
                with closing(request.urlopen(xml_url)) as xml_gz_file:
                    with open(xml_file, 'wb') as f:
                        with gzip.open(xml_gz_file, 'rb') as gzip_file:
                            shutil.copyfileobj(gzip_file, f)
                ElementTree.parse(xml_file)  # Just for validation purposes
            except (ConnectionResetError, ElementTree.ParseError) as e:
                print(f"Retrying download, reason: {e}")
                self.__download_xml(pdb_id)
            print(f"PDB to Uniprot XML downloaded successfully. Path {xml_file}")
        return xml_file

    def __parse_sifts(self, sifts_file):
        """
        Parse Sifts files and generate PdbChain objects
        """
        residues = {}  # { chain: { uniprot_id: { pdb_res_id : mapping_pos_index } } }
        nob = {}

        xml_entry = SiftsParser().parse(sifts_file)
        for ent in xml_entry.entities:
            for seg in ent.segments:
                for res in seg.residues:
                    if res.pdbresnum:
                        residues.setdefault(res.chainid, {}).setdefault(res.uniprotaccession, {})
                        residues[res.chainid][res.uniprotaccession][res.pdbresnum] = res.uniprotpos
                    if res.notobserved:
                        nob.setdefault(res.chainid, [])
                        nob[res.chainid].append(res.uniprotpos)
                    if res.uniprotaccession:
                        print("{},{},{},{},{},{},{},{},{},{},{},{},{}".format(
                            res.uniprotaccession, res.pdbid, res.chainid,
                            res.uniprotpos, res.pdbresnum, res.naturalpos, res.uniprotresname,
                            res.pdbresname, res.seqresname, res.notobserved, res.conflict,
                            res.modified, res.ss))

        return residues, nob


class SiftsEntry:
    def __init__(self, id):
        self.__id = id
        self.__entities = []

    @property
    def id(self):
        return self.__id

    @id.setter
    def id(self, id):
        self.__id = id

    @property
    def entities(self):
        return self.__entities

    @entities.setter
    def entities(self, entities):
        self.__entities = entities

    def addentity(self, entity):
        self.__entities.append(entity)

    def __getitem__(self, index):
        return self.__entities[index]

    def __setitem__(self, index, value):
        self.__entities[index] = value

    def __repr__(self):
        return "SiftsEntry {0}".format(self.__id)


class SiftsEntity:
    def __init__(self, id, type):
        self.__type = type
        self.__id = id
        self.__segments = []

    @property
    def type(self):
        return self.__type

    @type.setter
    def type(self, type):
        self.__type = type

    @property
    def id(self):
        return self.__id

    @id.setter
    def id(self, id):
        self.__id = id

    @property
    def segments(self):
        return self.__segments

    @segments.setter
    def segments(self, segments):
        self.__segments = segments

    def addsegment(self, segment):
        self.__segments.append(segment)

    def __getitem__(self, index):
        return self.__segments[index]

    def __setitem__(self, index, value):
        self.__segments[index] = value

    def __repr__(self):
        return "SiftsEntity (id: {0}, type:{1})".format(self.__id, self.__type)


class SiftsSegment:
    def __init__(self, id, start, end):
        self.__id = id
        self.__start = start
        self.__end = end
        self.__residues = []

    @property
    def id(self):
        return self.__id

    @id.setter
    def id(self, id):
        self.__id = id

    @property
    def start(self):
        return self.__start

    @start.setter
    def start(self, start):
        self.__start = start

    @property
    def end(self):
        return self.__end

    @end.setter
    def end(self, end):
        self.__end = end

    @property
    def residues(self):
        return self.__residues

    @residues.setter
    def residues(self, residues):
        self.__residues = residues

    def addresidue(self, residue):
        self.__residues.append(residue)

    def __getitem__(self, index):
        return self.__residues[index]

    def __setitem__(self, index, value):
        self.__residues[index] = value

    def __repr__(self):
        return "SiftsSegment (id: {0}, start:{1}, end: {2})".format(self.__id, self.__start, self.__end)


class SiftsResidue:
    def __init__(self):
        self.__pdbresnum = ""
        self.__pdbresname = ""
        self.__chainid = ""
        self.__uniprotresname = ""
        self.__uniprotpos = None
        self.__naturalpos = None
        self.__seqresname = ""
        self.__pdbid = ""
        self.__uniprotaccession = ""
        self.__ss = ""
        self.__notobserved = False
        self.__conflict = False
        self.__modified = False

    @property
    def pdbresnum(self):
        return self.__pdbresnum

    @pdbresnum.setter
    def pdbresnum(self, pdbresnum):
        self.__pdbresnum = pdbresnum

    @property
    def pdbresname(self):
        return self.__pdbresname

    @pdbresname.setter
    def pdbresname(self, pdbresname):
        self.__pdbresname = pdbresname

    @property
    def chainid(self):
        return self.__chainid

    @chainid.setter
    def chainid(self, chainid):
        self.__chainid = chainid

    @property
    def uniprotresname(self):
        return self.__uniprotresname

    @uniprotresname.setter
    def uniprotresname(self, uniprotresname):
        self.__uniprotresname = uniprotresname

    @property
    def uniprotpos(self):
        return self.__uniprotpos

    @uniprotpos.setter
    def uniprotpos(self, uniprotpos):
        self.__uniprotpos = uniprotpos

    @property
    def naturalpos(self):
        return self.__naturalpos

    @naturalpos.setter
    def naturalpos(self, naturalpos):
        self.__naturalpos = naturalpos

    @property
    def seqresname(self):
        return self.__seqresname

    @seqresname.setter
    def seqresname(self, seqresname):
        self.__seqresname = seqresname

    @property
    def pdbid(self):
        return self.__pdbid

    @pdbid.setter
    def pdbid(self, pdbid):
        self.__pdbid = pdbid

    @property
    def uniprotaccession(self):
        return self.__uniprotaccession

    @uniprotaccession.setter
    def uniprotaccession(self, uniprotaccession):
        self.__uniprotaccession = uniprotaccession

    @property
    def notobserved(self):
        return self.__notobserved

    @notobserved.setter
    def notobserved(self, notobserved):
        self.__notobserved = notobserved

    @property
    def ss(self):
        return self.__ss

    @ss.setter
    def ss(self, ss):
        self.__ss = ss

    @property
    def modified(self):
        return self.__modified

    @modified.setter
    def modified(self, modified):
        self.__modified = modified

    @property
    def conflict(self):
        return self.__conflict

    @conflict.setter
    def conflict(self, conflict):
        self.__conflict = conflict


class SiftsParsingError(Exception):
    pass


class SiftsParser:
    def __init__(self):
        self.__SIFTS_TAG_PREFIX = "{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}"

    def parse(self, filepath):
        try:
            f = open(filepath)
        except IOError:
            raise SiftsParsingError("{0} does not exist".format(filepath))
        root = ElementTree.parse(f).getroot()
        pdbid = root.attrib["dbAccessionId"]

        sen = SiftsEntry(pdbid)

        for entity in root.findall('{0}entity'.format(self.__SIFTS_TAG_PREFIX)):

            type = entity.attrib["type"]
            id = entity.attrib["entityId"]

            se = SiftsEntity(type, id)

            sen.addentity(se)

            for segment in entity.findall("{0}segment".format(self.__SIFTS_TAG_PREFIX)):

                id = segment.attrib["segId"]
                start = segment.attrib["start"]
                end = segment.attrib["end"]

                ss = SiftsSegment(id, start, end)

                se.addsegment(ss)

                for residue in segment.findall("{0}listResidue/{0}residue".format(self.__SIFTS_TAG_PREFIX)):

                    sr = SiftsResidue()

                    ss.addresidue(sr)

                    sr.naturalpos = int(residue.attrib["dbResNum"])
                    sr.seqresname = residue.attrib["dbResName"]

                    for detail in residue.findall("{0}residueDetail".format(self.__SIFTS_TAG_PREFIX)):

                        if detail.attrib["property"] == "codeSecondaryStructure":
                            sr.ss = detail.text
                        elif detail.attrib["property"] == "Annotation":
                            ann = detail.text.strip().replace("\n          ", "")
                            if ann == "Not_Observed":
                                sr.notobserved = True
                            elif ann == "Conflict":
                                sr.conflict = True
                            elif ann == "PDB modified":
                                sr.modified = True

                    links = residue.findall("{0}crossRefDb".format(self.__SIFTS_TAG_PREFIX))

                    for link in links:
                        dbsrc = link.attrib["dbSource"]
                        dbacc = link.attrib["dbAccessionId"]
                        dbresnum = link.attrib["dbResNum"]
                        dbresname = link.attrib["dbResName"]

                        if dbsrc == "PDB":
                            dbchainid = link.attrib["dbChainId"]

                            sr.pdbid = dbacc
                            sr.chainid = dbchainid
                            sr.pdbresnum = dbresnum
                            sr.pdbresname = dbresname
                        elif dbsrc == "UniProt":
                            sr.uniprotaccession = dbacc
                            sr.uniprotpos = int(dbresnum)
                            sr.uniprotresname = dbresname

        return sen
