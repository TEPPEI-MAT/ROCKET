import math
from .dna import DNA_ENERGIES
from typing import List, Tuple, Union
from .fold import dg
import argparse
from Bio import SeqIO
import os.path

parser = argparse.ArgumentParser(description="Rational Oligonucleotide design Calculated with Kinetic parameter for Enhanced in vitro Transcription")

parser.add_argument(
    "-s", "--sequence", type=str, help="tRNA sequence to transcribe"
)
parser.add_argument(
    "-d", "--dir_path", type=str, help="Path to fasta file (Output as 'output.txt' in the same directory)"
)
parser.add_argument(
    "-T", "--maximum_tm", type=int, default=72, help="maximum Tm (default is 72 degrees C)"
)
parser.add_argument(
    "-t", "--minimum_tm", type=int, default=68,  help="minimum Tm (default is 68 degrees C)"
)
parser.add_argument(
    "-M", "--maximum_primer", type=int, default=100, help="maximum length of primer (default is 100 nt)"
)
parser.add_argument(
    "-m", "--minimum_primer", type=int, default=7,  help="minimum length of primer (default is 7 nt)"
)
parser.add_argument(
    "-p", "--precursor", action='store_true', help="Add precursor to the 5' end"
)
parser.add_argument(
    "-g", "--g_addition", action='store_true', help="Add a guanine base to the 5' end"
)
parser.add_argument(
    "-c", "--cca_addition", action='store_true', help="Add CCA sequence to the 3' end"
)

args = parser.parse_args()

"""
The ROCKET is wrote to use on python or IDE such as VScode, jupyternote and Pycharm.
"""

def ROCKET(seq:str, g=False, c=False, p=False, min_nt:int=7, max_nt:int=100, min_tm:float = 68, max_tm:float = 72):
    if g == True:
        seq = "G" + seq
    else:
        pass
    if c == True:
        seq = seq + "CCA"
    else:
        pass

    def RNA_Cuter(trancated_seq,start_nt: int = 1,end_nt: int = len(seq)) -> List:
        for i in range(start_nt - 1, end_nt):
            NN = [seq[i:n+i+1] for n in range(0, len(seq)-i, 1)]
            trancated_seq.extend(NN)
        return trancated_seq

    def Tm_Calculator(seq:str, ct:float = 100, Na:float = 1) -> Union[float, str] :
        seq = seq.upper()
        seq = seq.replace("U", "T")
        NN = [seq[n:n+2] for n in range(len(seq)-1)]
        seq_dH = [DNA_ENERGIES.NN2[n][0] for n in NN]
        seq_dS = [DNA_ENERGIES.NN2[n][1] for n in NN]

        if len(NN) != 0:
            if NN[0][0] == "G" or NN[0][0] == "C":
                S_sum = sum(seq_dS) - 16.8

                if NN[-1][1] == "A" or NN[-1][1] ==  "T":
                    H_sum = sum(seq_dH) + 2.2
                    S_sum = S_sum + 6.9

                else:
                    H_sum = sum(seq_dH)

            elif NN[0][0] == "A" or NN[0][0] ==  "T":
                S_sum = sum(seq_dS) - 20.1

                if NN[-1][1] == "A" or NN[-1][1] == "T":
                    H_sum = sum(seq_dH) + 2.2
                    S_sum = S_sum + 6.9

                else:
                    H_sum = sum(seq_dH)

            else:
                S_sum = sum(seq_dS)

                if NN[-1][1] == "A" or NN[-1][1] == "T":
                    H_sum = sum(seq_dH) + 2.2
                    S_sum = S_sum + 6.9

                else:
                    H_sum = sum(seq_dH)
        else:
            H_sum = sum(seq_dH)
            S_sum = sum(seq_dS)

        R = 8.314462618
        Ct = ct * 1e-6
        if Na == 1 :
            Tm = H_sum * (10 ** 3) / (S_sum + R * math.log(Ct * 2 / 4)) - 273

        elif 0.05 <= Na <= 1.1 :
            S_sum = S_sum + (0.368 * (len(seq) - 1)) * math.log(Na)
            Tm = H_sum * (10 ** 3) / (S_sum + R * math.log(Ct * 2 / 4)) - 273

        else:
            raise ValueError("Accuracy is guaranteed only during 0.05 <= Na+ <= 1.1")

        return Tm

    def Index_Get_Tm(lst:List[float]) -> List:
        return [i for i , _x in enumerate(lst) if min_tm <= _x <= max_tm]

    def Index_Get_GC(lst:List[str]) -> List:
        return [i for i , _x in enumerate(lst) if (_x[0] == "G" or _x[0] == "C") and (_x[-1] == "G" or _x[-1] == "C")]

    def Reverse_Complement_Strand(seq:str) -> str:
        translation = seq.maketrans("AUTGC", "TAACG")
        template = seq[::-1]
        complement = template.translate(translation)
        return complement

    def T7_Promoter_Attach(seq:str, aneeling_nt:str) -> Tuple[str, str]:
        seq = seq.upper()
        seq = seq.replace("U","T")
        anl_point = seq.find(aneeling_nt.replace("U","T"))
        if p == True:
            Forward_Primer = "GCCTAATACGACTCACTATAGGGAGACCACAACGGTTTCCCTCTAGA" + seq[:(anl_point + len(aneeling_nt))]
        else:
            Forward_Primer = "GCCTAATACGACTCACTATA" + seq[:(anl_point + len(aneeling_nt))]
        Reverse_Primer = Reverse_Complement_Strand(seq[anl_point:])
        return Forward_Primer, Reverse_Primer

    def ANNEALING_AREA() -> List[str]:

        trancated_seq = []

        Trancated_seq = RNA_Cuter(trancated_seq)
        Tm = [Tm_Calculator(i) for i in Trancated_seq]

        narrowed_seq = [Trancated_seq[i] for i in Index_Get_Tm(Tm)]

        High_narrowed_seq = [narrowed_seq[i] for i in Index_Get_GC(narrowed_seq)]
        return High_narrowed_seq

    High_narrowed_seq = ANNEALING_AREA()

    primers = [(T7_Promoter_Attach(seq, High_narrowed_seq[i])) for i in range(len(High_narrowed_seq))]
    Primers = []
    for i in range(len(primers)):
        if min_nt <= len((primers[i][0])) <= max_nt and min_nt <= len(primers[i][1]) <= max_nt:
            Primers.append(primers[i])

    dG = [dg(Primers[i][0]) + dg(Primers[i][1]) for i in range(len(Primers))]
    top = [x for x, _x in enumerate(dG) if _x == max(dG)]
    l = []
    if top == 1:
        print(">" + "Forward")
        print(Primers[top[0][0]])
        print(">" + "Reverse")
        print(Primers[top[0][1]])
        #print("Total dG: " + str(max(dG)) + "kcal/mol")

    else:
        for i in range(len(top)):
            a = len(Primers[top[i]][0]) + len(Primers[top[i]][1])
            l.append(a)

        b = [i for i, x in enumerate(l) if x == min(l)]
        if len(b) == 1:
            print(">" + "Forward")
            print(Primers[top[b[0]]][0])
            print(">" + "Reverse")
            print(Primers[top[b[0]]][1])
            #print("Total dG: " + str(max(dG)) + "kcal/mol")

        else:
            for i in range(len(b)):
                c = abs(len(Primers[top[i]][0]) - len(Primers[top[i]][1]))
                l.append(c)
            d = [i for i, x in enumerate(l) if x == min(l)]
            if len(d) == 1:
                print(">" + "Forward")
                print(Primers[top[b[0]]][0])
                print(">" + "Reverse")
                print(Primers[top[b[0]]][1])
                #print("Total dG: " + str(max(dG)) + "kcal/mol")

            else:
                for i in range(len(d)):
                    print(">" + "Forward")
                    print(Primers[top[0][0]])
                    print(">" + "Reverse")
                    print(Primers[top[0][1]])
                    #print("Total dG: " + str(max(dG)) + "kcal/mol")

def run():
    if args.dir_path == None:
        seq = args.sequence
        if args.g_addition:
            seq = "G" + seq
        else:
            pass
        if args.cca_addition:
            seq = seq + "CCA"
        else:
            pass


        def RNA_Cuter(trancated_seq,start_nt: int = 1,end_nt: int = len(seq)) -> List:
            for i in range(start_nt - 1, end_nt):
                NN = [seq[i:n+i+1] for n in range(0, len(seq)-i, 1)]
                trancated_seq.extend(NN)
            return trancated_seq

        def Tm_Calculator(seq:str, ct:float = 100, Na:float = 1) -> Union[float, str] :
            seq = seq.upper()
            seq = seq.replace("U", "T")
            NN = [seq[n:n+2] for n in range(len(seq)-1)]
            seq_dH = [DNA_ENERGIES.NN2[n][0] for n in NN]
            seq_dS = [DNA_ENERGIES.NN2[n][1] for n in NN]

            if len(NN) != 0:
                if NN[0][0] == "G" or NN[0][0] == "C":
                    S_sum = sum(seq_dS) - 16.8

                    if NN[-1][1] == "A" or NN[-1][1] ==  "T":
                        H_sum = sum(seq_dH) + 2.2
                        S_sum = S_sum + 6.9

                    else:
                        H_sum = sum(seq_dH)

                elif NN[0][0] == "A" or NN[0][0] ==  "T":
                    S_sum = sum(seq_dS) - 20.1

                    if NN[-1][1] == "A" or NN[-1][1] == "T":
                        H_sum = sum(seq_dH) + 2.2
                        S_sum = S_sum + 6.9

                    else:
                        H_sum = sum(seq_dH)

                else:
                    S_sum = sum(seq_dS)

                    if NN[-1][1] == "A" or NN[-1][1] == "T":
                        H_sum = sum(seq_dH) + 2.2
                        S_sum = S_sum + 6.9

                    else:
                        H_sum = sum(seq_dH)
            else:
                H_sum = sum(seq_dH)
                S_sum = sum(seq_dS)

            R = 8.314462618
            Ct = ct * 1e-6
            if Na == 1 :
                Tm = H_sum * (10 ** 3) / (S_sum + R * math.log(Ct * 2 / 4)) - 273

            elif 0.05 <= Na <= 1.1 :
                S_sum = S_sum + (0.368 * (len(seq) - 1)) * math.log(Na)
                Tm = H_sum * (10 ** 3) / (S_sum + R * math.log(Ct * 2 / 4)) - 273

            else:
                raise ValueError("Accuracy is guaranteed only during 0.05 <= Na+ <= 1.1")

            return Tm

        def Index_Get_Tm(lst:List[float]) -> List:
            return [i for i , _x in enumerate(lst) if args.minimum_tm <= _x <= args.maximum_tm]

        def Index_Get_GC(lst:List[str]) -> List:
            return [i for i , _x in enumerate(lst) if (_x[0] == "G" or _x[0] == "C") and (_x[-1] == "G" or _x[-1] == "C")]

        def Reverse_Complement_Strand(seq:str) -> str:
            translation = seq.maketrans("AUTGC", "TAACG")
            template = seq[::-1]
            complement = template.translate(translation)
            return complement

        def T7_Promoter_Attach(seq:str, aneeling_nt:str) -> Tuple[str, str]:
            seq = seq.upper()
            seq = seq.replace("U","T")
            anl_point = seq.find(aneeling_nt.replace("U","T"))
            if args.precursor:
                Forward_Primer = "GCCTAATACGACTCACTATAGGGAGACCACAACGGTTTCCCTCTAGA" + seq[:(anl_point + len(aneeling_nt))]
            else:
                Forward_Primer = "GCCTAATACGACTCACTATA" + seq[:(anl_point + len(aneeling_nt))]
            Reverse_Primer = Reverse_Complement_Strand(seq[anl_point:])
            return Forward_Primer, Reverse_Primer

        def ANNEALING_AREA() -> List[str]:

            trancated_seq = []

            Trancated_seq = RNA_Cuter(trancated_seq)
            Tm = [Tm_Calculator(i) for i in Trancated_seq]

            narrowed_seq = [Trancated_seq[i] for i in Index_Get_Tm(Tm)]

            High_narrowed_seq = [narrowed_seq[i] for i in Index_Get_GC(narrowed_seq)]
            return High_narrowed_seq

        High_narrowed_seq = ANNEALING_AREA()

        primers = [(T7_Promoter_Attach(seq, High_narrowed_seq[i])) for i in range(len(High_narrowed_seq))]
        Primers = []
        for i in range(len(primers)):
            if args.minimum_primer <= len((primers[i][0])) <= args.maximum_primer and args.minimum_primer <= len(primers[i][1]) <= args.maximum_primer:
                Primers.append(primers[i])

        dG = [dg(Primers[i][0]) + dg(Primers[i][1]) for i in range(len(Primers))]
        top = [x for x, _x in enumerate(dG) if _x == max(dG)]
        l = []
        la = []
        if top == 1:
            print(">" + "Forward")
            print(Primers[top[0][0]])
            print(">" + "Reverse")
            print(Primers[top[0][1]])
            #print("Total dG: " + str(max(dG)) + "kcal/mol")

        else:
            for i in range(len(top)):
                a = len(Primers[top[i]][0]) + len(Primers[top[i]][1])
                l.append(a)

            b = [i for i, x in enumerate(l) if x == min(l)]
            if len(b) == 1:
                print(">" + "Forward")
                print(Primers[top[b[0]]][0])
                print(">" + "Reverse")
                print(Primers[top[b[0]]][1])
                #print("Total dG: " + str(max(dG)) + "kcal/mol")

            else:
                for i in range(len(b)):
                    c = abs(len(Primers[top[i]][0]) - len(Primers[top[i]][1]))
                    l.append(c)
                d = [i for i, x in enumerate(l) if x == min(l)]
                if len(d) == 1:
                    print(">" + "Forward")
                    print(Primers[top[b[0]]][0])
                    print(">" + "Reverse")
                    print(Primers[top[b[0]]][1])
                    #print("Total dG: " + str(max(dG)) + "kcal/mol")

                else:
                    print(">" + "Forward")
                    print(Primers[top[b[0]]][0])
                    print(">" + "Reverse")
                    print(Primers[top[b[0]]][1])
                    #print("Total dG: " + str(max(dG)) + "kcal/mol")

    else:
        Seq = []
        iD = []


        def read_fasta_file(filename):
            sequences = []

            for record in SeqIO.parse(filename, "fasta"):
                sequences.append(record)

            return sequences


        sequence = read_fasta_file(args.dir_path)

        for i in sequence:
            s = str(i.seq)
            Seq.append(s)
            o = str(i.id)
            iD.append(o)

        for T in range(len(Seq)):
            seq = Seq[T]
            if args.g_addition:
                seq = "G" + seq
            else:
                pass
            if args.cca_addition:
                seq = seq + "CCA"
            else:
                pass


            def RNA_Cuter(trancated_seq, start_nt: int = 1, end_nt: int = len(seq)) -> List:
                for i in range(start_nt - 1, end_nt):
                    NN = [seq[i:n + i + 1] for n in range(0, len(seq) - i, 1)]
                    trancated_seq.extend(NN)
                return trancated_seq


            def Tm_Calculator(seq: str, ct: float = 100, Na: float = 1) -> Union[float, str]:
                seq = seq.upper()
                seq = seq.replace("U", "T")
                NN = [seq[n:n + 2] for n in range(len(seq) - 1)]
                seq_dH = [DNA_ENERGIES.NN2[n][0] for n in NN]
                seq_dS = [DNA_ENERGIES.NN2[n][1] for n in NN]

                if len(NN) != 0:
                    if NN[0][0] == "G" or NN[0][0] == "C":
                        S_sum = sum(seq_dS) - 16.8

                        if NN[-1][1] == "A" or NN[-1][1] == "T":
                            H_sum = sum(seq_dH) + 2.2
                            S_sum = S_sum + 6.9

                        else:
                            H_sum = sum(seq_dH)

                    elif NN[0][0] == "A" or NN[0][0] == "T":
                        S_sum = sum(seq_dS) - 20.1

                        if NN[-1][1] == "A" or NN[-1][1] == "T":
                            H_sum = sum(seq_dH) + 2.2
                            S_sum = S_sum + 6.9

                        else:
                            H_sum = sum(seq_dH)

                    else:
                        S_sum = sum(seq_dS)

                        if NN[-1][1] == "A" or NN[-1][1] == "T":
                            H_sum = sum(seq_dH) + 2.2
                            S_sum = S_sum + 6.9

                        else:
                            H_sum = sum(seq_dH)
                else:
                    H_sum = sum(seq_dH)
                    S_sum = sum(seq_dS)

                R = 8.314462618
                Ct = ct * 1e-6
                if Na == 1:
                    Tm = H_sum * (10 ** 3) / (S_sum + R * math.log(Ct * 2 / 4)) - 273

                elif 0.05 <= Na <= 1.1:
                    S_sum = S_sum + (0.368 * (len(seq) - 1)) * math.log(Na)
                    Tm = H_sum * (10 ** 3) / (S_sum + R * math.log(Ct * 2 / 4)) - 273

                else:
                    raise ValueError("Accuracy is guaranteed only during 0.05 <= Na+ <= 1.1")

                return Tm


            def Index_Get_Tm(lst: List[float]) -> List:
                return [i for i, _x in enumerate(lst) if args.minimum_tm <= _x <= args.maximum_tm]


            def Index_Get_GC(lst: List[str]) -> List:
                return [i for i, _x in enumerate(lst) if
                        (_x[0] == "G" or _x[0] == "C") and (_x[-1] == "G" or _x[-1] == "C")]


            def Reverse_Complement_Strand(seq: str) -> str:
                translation = seq.maketrans("AUTGC", "TAACG")
                template = seq[::-1]
                complement = template.translate(translation)
                return complement


            def T7_Promoter_Attach(seq: str, aneeling_nt: str) -> Tuple[str, str]:
                seq = seq.upper()
                seq = seq.replace("U", "T")
                anl_point = seq.find(aneeling_nt.replace("U", "T"))
                if args.precursor:
                    Forward_Primer = "GCCTAATACGACTCACTATAGGGAGACCACAACGGTTTCCCTCTAGA" + seq[
                                                                                         :(anl_point + len(aneeling_nt))]
                else:
                    Forward_Primer = "GCCTAATACGACTCACTATA" + seq[:(anl_point + len(aneeling_nt))]
                Reverse_Primer = Reverse_Complement_Strand(seq[anl_point:])
                return Forward_Primer, Reverse_Primer


            def ANNEALING_AREA() -> List[str]:

                trancated_seq = []

                Trancated_seq = RNA_Cuter(trancated_seq)
                Tm = [Tm_Calculator(i) for i in Trancated_seq]

                narrowed_seq = [Trancated_seq[i] for i in Index_Get_Tm(Tm)]

                High_narrowed_seq = [narrowed_seq[i] for i in Index_Get_GC(narrowed_seq)]
                return High_narrowed_seq


            High_narrowed_seq = ANNEALING_AREA()

            primers = [(T7_Promoter_Attach(seq, High_narrowed_seq[i])) for i in range(len(High_narrowed_seq))]
            Primers = []
            for i in range(len(primers)):
                if args.minimum_primer <= len((primers[i][0])) <= args.maximum_primer and args.minimum_primer <= len(primers[i][1]) <= args.maximum_primer:
                    Primers.append(primers[i])

            dG = [dg(Primers[i][0]) + dg(Primers[i][1]) for i in range(len(Primers))]
            top = [x for x, _x in enumerate(dG) if _x == max(dG)]
            l = []
            la = []
            if top == 1:
                path = os.path.dirname(args.dir_path)
                with open(path + "\output.txt", "a") as f:
                    f.write(">" + iD[T] + "Forward\n")
                    f.write(Primers[top[0][0]] + "\n")
                    f.write(">" + iD[T] + "Reverse\n")
                    f.write(Primers[top[0][1]] + "\n")

            else:
                for i in range(len(top)):
                    a = len(Primers[top[i]][0]) + len(Primers[top[i]][1])
                    l.append(a)

                b = [i for i, x in enumerate(l) if x == min(l)]
                if len(b) == 1:
                    path = os.path.dirname(args.dir_path)
                    with open(path + "\output.txt", "a") as f:
                        f.write(">" + iD[T] + "_Forward\n")
                        f.write(Primers[top[b[0]]][0] + "\n")
                        f.write(">" + iD[T] + "_Reverse\n")
                        f.write(Primers[top[b[0]]][1] + "\n")

                else:
                    for i in range(len(b)):
                        c = abs(len(Primers[top[i]][0]) - len(Primers[top[i]][1]))
                        l.append(c)
                    d = [i for i, x in enumerate(l) if x == min(l)]
                    if len(d) == 1:
                        path = os.path.dirname(args.dir_path)
                        with open(path + "\output.txt", "a") as f:
                            f.write(">" + iD[T] + "_Forward\n")
                            f.write(Primers[top[b[0]]][0] + "\n")
                            f.write(">" + iD[T] + "_Reverse\n")
                            f.write(Primers[top[b[0]]][1] + "\n")

                    else:
                        for i in range(len(d)):
                            print(iD[T])
                            print(Primers[top[i]])
                            print("Total dG: " + str(max(dG)) + "kcal/mol")

if __name__ == "__main__":
    run()
    print("completed")
