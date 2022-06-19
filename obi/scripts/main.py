def main():
    print(f'There are the available commands, each one has its own parameters. Use obi <command> --help to learn more '
          f'about the usage\n')
    print(f'obi fetch-db                  Downloads Swissprot DB and Uniprot <-> PDB mappings from NCBI FTP server')
    print(f'obi run [<args>]              Runs the whole analysis or first part depending on the arguments')
    print(f'obi analysis [<args>]         Runs positive selection analysis, which includes running Hyphy and getting '
          f'final results')
    print(f'obi resume-analysis [<args>]  Runs the whole analysis or first part depending on the arguments')


if __name__ == "__main__":
    main()
