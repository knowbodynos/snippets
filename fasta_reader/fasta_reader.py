import os

class FastaReader:
    def __init__(self, fasta, nb_line=60):
        self._chrm = chrm
        self._start = start
        self._end = end
        self._nb_line = nb_line
        self._fasta_stream, self._bytes_like = self._get_fasta_stream(fasta)

    def _get_fasta_stream(fasta):
        if isinstance(fasta, str):
            if os.path.isfile(fasta):
                with open(fasta, 'r') as fasta_stream:
                    char = fasta_stream.read(1)
                    bytes_like = hasattr(char, decode)
                    fasta_stream.seek(0)
                    return fasta_stream, False
            else:
                raise ValueError("{} not a valid file path.".format(fasta))
        elif hasattr(fasta, 'read') and hasattr(fasta, 'readline') and \
             hasattr(fasta, 'seek') and hasattr(fasta, 'tell'):
            char = fasta.read(1)
            bytes_like = hasattr(char, decode)
            fasta.seek(0)
            return fasta, bytes_like
        else:
            raise ValueError("You must provide a valid file path or a file stream.".format(fasta))

    def set_fasta(fasta):
        self._fasta_stream, self._bytes_like = self._get_fasta_stream(fasta)

    def set_nb_line(nb_line):
        self._nb_line = nb_line

    def seq(chrm, start, end, nb_line=None):
        if nb_line is None:
            nb_line = self._nb_line
        if self._bytes_like:
            sep = b':'
            eol = b'\n'
            seq = b''
            chrm = chrm.encode()
        else:
            sep = ':'
            eol = '\n'
            seq = ''
        pos_bkp = self._fasta_stream.tell()
        self._fasta_stream.seek(0)
        t_chrm = None
        # Loop until reaching specified chromosome
        while t_chrm != chrm:
            if not t_chrm is None:
                # Seek until next chromosome
                ## t_nb // nb_line is the number of newline characters, 
                ## add last newline before next chromosome header
                self._fasta_stream.seek(self._fasta_stream.tell() + t_nb + (t_nb // nb_line) + 1)
            line = self._fasta_stream.readline()
            # Grab chromosome name and number of bases from chromosome header line
            t_chrm, _, t_nb = line.split(sep)[3:6]
            t_nb = int(t_nb)
        i = start - 1
        file_start = self._fasta_stream.tell() + i + (i // nb_line)
        # Seek to just before specified sequence 
        self._fasta_stream.seek(file_start)
        # Loop until end of specified sequence
        while i < end:
            # If end of sequence is reached on current line
            if i // nb_line == (end - 1) // nb_line:
                line = self._fasta_stream.read(end - i)
            else:
                line = self._fasta_stream.readline()
            line = line.rstrip(eol)
            seq += line
            i += len(line)
        self._fasta_stream.seek(pos_bkp)
        if self._bytes_like:
            return(seq.decode())
        else:
            return(seq)