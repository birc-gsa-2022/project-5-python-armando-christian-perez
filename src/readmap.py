import argparse
import sys
from dataclasses import dataclass,field

def main():
    argparser = argparse.ArgumentParser(
        description="Readmapper",
        usage="\n\treadmap -p genome\n\treadmap -d dist genome reads"
    )
    argparser.add_argument(
        "-p", action="store_true",
        help="preprocess the genome."
    )
    argparser.add_argument(
        "-d", type=int, metavar="integer",
        default=1, help="max edit distance."
    )
    argparser.add_argument(
        "genome",
        help="Simple-FASTA file containing the genome.",
        type=argparse.FileType('r')
    )
    argparser.add_argument(
        "reads", nargs="?",
        help="Simple-FASTQ file containing the reads.",
        type=argparse.FileType('r')
    )
    args = argparser.parse_args()
    name = args.genome.name.replace(".fa", "")
    if args.p:
        create_BWT_structures(args.genome, name)
    else:
        fasta_dict, fastq_dict = fasta_translator(args.genome), fastq_translator(args.reads)
        SAM = matches_to_SAM(fasta_dict, fastq_dict, name, args.d)
        print_SAM(SAM)
        if args.reads is None:
            argparser.print_help()
            sys.exit(1)

@dataclass
class node:
    # A node in a suffix tree
    index: list = field(default_factory=lambda: [0, None])
    suffix_link: int = None
    parent: int = None # Batman
    children: dict = field(default_factory = dict)
    self_pointer: int = 0
    leaf_start: int = None

def fasta_translator(input_file):
    output_dict = {}
    start = True
    for i in input_file:
        i = i.strip()
        if len(i) == 0:
            continue
        if i[0] == ">":
            if  not start:
                output_dict[name] = seq + "$"
            name = i[1:].strip()
            seq = ""
            if start:
                start = False
        elif start:
            continue
        else:
            seq += i
    output_dict[name] = seq + "$"
    return output_dict

def fastq_translator(input_file):
    output_dict = {}
    for i in input_file:
        if i[0] == "@":
            name = i[1:].strip()
        else:
            seq = i.strip()
            output_dict[name] = seq
    return(output_dict)

def suffix_tree(string):
    remaining = 0
    active_edge = None
    edge_length = 0
    tree_list = [node(index = None)] #create root node
    active_node = tree_list[0]
    string = string
    for i in range(len(string)):
        remaining += 1
        current_char = string[i]
        suffix_link_update = None
        while remaining > 0: # Not sure if this is spaghetti code, if there are some good code practices i break here please lemme know
            if edge_length: # there is an active edge
                edge_char = string[active_edge]
                matched_node = tree_list[active_node.children[edge_char]]

                if matched_node.index[1] is not None and matched_node.index[1] < (matched_node.index[0] + edge_length): # We have to jump an internal node
                    internal_node_length = matched_node.index[1] - matched_node.index[0] + 1

                    if internal_node_length < edge_length: # We have to go to a child of the internal node to begin matching
                        active_node = matched_node
                        active_edge += internal_node_length
                        edge_length -= internal_node_length
                        continue # we have the code to care of the rest elsewhere, no need to write it twice in a while loop

                    else: # We can begin matching with the start of the internal nodes edges
                        if current_char in matched_node.children.keys(): # we can extend from the internal node
                            active_node = matched_node
                            matched_node = tree_list[active_node.children[current_char]]
                            active_edge = matched_node.index[0]
                            edge_length = 1
                            if suffix_link_update:
                                tree_list[suffix_link_update].suffix_link = active_node.self_pointer
                            break

                        else: # we can't extend from the internal node

                            newnode = node(index = [i, None], parent = matched_node.self_pointer, self_pointer = len(tree_list),
                                leaf_start = i - remaining + 1)
                            tree_list.append(newnode)
                            tree_list[matched_node.self_pointer].children[current_char] = len(tree_list) - 1
                            if suffix_link_update:
                                tree_list[suffix_link_update].suffix_link = matched_node.self_pointer
                            suffix_link_update = matched_node.self_pointer
                            remaining -= 1
                            if active_node.suffix_link is not None:
                                active_node = tree_list[active_node.suffix_link]

                            else:
                                active_edge += 1
                                edge_length -= 1
                                active_node == tree_list[0]

                else: # We don't have to jump an internal node
                    match = string[matched_node.index[0] + edge_length] # the character we are matching against
                    if current_char == match: # there is a match
                        edge_length += 1
                        break

                    else: # there is no match
                        newnode = node(index = [matched_node.index[0], matched_node.index[0] + edge_length - 1], parent = active_node.self_pointer,
                            self_pointer = len(tree_list), suffix_link = 0)
                            # creates a new internal node which has the old matched node or leaf as its child
                        internal_node_index = newnode.self_pointer
                        newnode.children[string[matched_node.index[0] + edge_length]] = matched_node.self_pointer
                        tree_list[matched_node.self_pointer].parent = internal_node_index

                        tree_list.append(newnode)
                        tree_list[active_node.self_pointer].children[string[newnode.index[0]]] = internal_node_index
                        tree_list[matched_node.self_pointer].index[0] = newnode.index[1] + 1

                        if suffix_link_update:
                            tree_list[suffix_link_update].suffix_link = internal_node_index
                        suffix_link_update = internal_node_index

                        newnode = node(index = [i, None], parent = internal_node_index, self_pointer = len(tree_list),
                        leaf_start = i - remaining + 1) # and creates the new leaf
                        tree_list.append(newnode)
                        tree_list[internal_node_index].children[current_char] = newnode.self_pointer
                        remaining -= 1
                        if active_node.suffix_link is not None:
                            active_node = tree_list[active_node.suffix_link]
                        else:
                            active_edge += 1
                            edge_length -= 1
                            active_node == tree_list[0]
            else: #There is no active edge
                
                assert active_node.self_pointer == 0
                if current_char in active_node.children.keys(): # We can extend from root
                    active_edge = tree_list[active_node.children[current_char]].index[0]
                    edge_length = 1
                    if suffix_link_update:
                        tree_list[suffix_link_update].suffix_link = active_node.self_pointer
                    break
                else: # we cant extend from root
                    newnode = node(index = [i, None], parent = active_node.self_pointer, self_pointer = len(tree_list), leaf_start = i - remaining + 1)
                    tree_list.append(newnode)
                    tree_list[active_node.self_pointer].children[current_char] = len(tree_list) - 1
                    remaining -= 1
                    if suffix_link_update:
                        tree_list[suffix_link_update].suffix_link = active_node.self_pointer
    return tree_list

def tree_to_SA(tree_list, string):
    stack = []
    SA = [None for i in string]
    count = 0
    stack.append(tree_list[0])
    while stack:
        if stack[-1].leaf_start is not None:
            node = stack.pop()
            SA[count] = node.leaf_start
            count += 1
        elif stack[-1].children:
            keys = [i for i in stack[-1].children.keys()]
            keys.sort()
            index = stack[-1].children.pop(keys[0])
            stack.append(tree_list[index])
        else:
            stack.pop()
    return SA

def SA_st(string):
    tree_list = suffix_tree(string)
    SA = tree_to_SA(tree_list, string)
    return SA

def compute_buckets(SA, reference):
    bucket_len = {}
    for i in SA:
        if reference[i - 1] in bucket_len:
            bucket_len[reference[i - 1]] += 1
        else:
            bucket_len[reference[i - 1]] = 1
    bucket_len = dict(sorted(bucket_len.items()))
    buckets = {}
    len_traversed = 0

    for key in bucket_len:
        temp = len_traversed
        len_traversed += bucket_len[key]
        buckets[key] = temp

    return buckets

def create_tracker(buckets, SA):
    tracker = {}
    for i in buckets:
        tracker[i] = [None for i in range(len(SA) + 1)]
    return tracker

def create_ranks(buckets, SA, reference):
    tracker = create_tracker(buckets, SA)
    for char in tracker:
        tracker[char][0] = 0
    for i in range(1, len(SA) + 1):
        for char in tracker:
            tracker[char][i] = tracker[char][i - 1]
        prev_bwt_char = reference[SA[i - 1] - 1]
        tracker[prev_bwt_char][i] += 1
    return tracker

def create_BWT_structures(ref_file, name):
    fasta_dict = fasta_translator(ref_file)

    SA_file = open(name + r"_SA.txt", "w")
    bucket_file = open(name + r"_buckets.txt", "w")
    rank_file = open(name + r"_rank.txt", "w")

    rev_SA_file = open(name + r"_rev_SA.txt", "w")
    rev_rank_file = open(name + r"_rev_rank.txt", "w")
    for key in fasta_dict:
        rev = fasta_dict[key][-2::-1] + "$"

        SA = SA_st(fasta_dict[key])
        SA_print = ",".join(str(i) for i in SA)
        SA_file.write(SA_print)
        
        SA_rev = SA_st(rev)
        SA_print_rev = ",".join(str(i) for i in SA_rev)
        rev_SA_file.write(SA_print)

        buckets = compute_buckets(SA, fasta_dict[key])
        for bucket in buckets:
            bucket_file.write(bucket + " " + str(buckets[bucket]) + "\t")

        rank = create_ranks(buckets, SA, fasta_dict[key])
        rank_chars = list(rank.keys())
        rank_file.write(",".join(rank_chars))
        for i in range(len(SA) + 1):
            rank_file.write("\t")
            row = [str(rank[char][i]) for char in rank_chars]
            rank_file.write(",".join(row))

        rev_rank = create_ranks(buckets, SA_rev, rev)
        rev_rank_chars = list(rank.keys())
        rev_rank_file.write(",".join(rev_rank_chars))
        for i in range(len(SA_rev) + 1):
            rev_rank_file.write("\t")
            row = [str(rev_rank[char][i]) for char in rev_rank_chars]
            rev_rank_file.write(",".join(row))

        rank_file.write("\n")
        SA_file.write("\n")
        bucket_file.write("\n")

        rev_rank_file.write("\n")
        rev_SA_file.write("\n")

    SA_file.close()
    bucket_file.close()
    rank_file.close()

    rev_SA_file.close()
    rev_rank_file.close()

def read_SA(line):
    SA = line.strip().split(",")
    SA = [int(sa) for sa in SA]
    return SA

def read_bucket(line):
    buckets = {}
    ite = line.strip().split("\t")
    for i in ite:
        info = i.split(" ")
        buckets[info[0]] = int(info[1])
    return buckets

def read_rank(line):
    ite = line.strip().split("\t")
    keys = ite[0].split(",")
    ranks = {key: [None for i in range(len(ite) - 1)] for key in keys}
    for i in range(1, len(ite)):
        row = [int(c) for c in ite[i].split(",")]
        for j, key in enumerate(keys):
            ranks[key][i - 1] = row[j]
    return ranks

def compute_D_table(read, rev_SA, buckets, rev_ranks):
    D_table = [None for i in read]
    upper, lower = 0, len(rev_SA)
    min_edits = 0

    for i, char in enumerate(read):
        upper = buckets[char] + rev_ranks[char][upper] # Will break if read char is not in genome
        lower = buckets[char] + rev_ranks[char][lower]
        if upper == lower:
            min_edits += 1
            upper, lower = 0, len(rev_SA)
        D_table[i] = min_edits
    return D_table

def Smoke_cigar(CIGAR):
    RTEC = ""
    prev = CIGAR[0]
    count = 1
    for char in CIGAR[1:]:
        if char == prev:
            count += 1
        elif count > 1:
            RTEC += (str(count) + prev)
            count = 1
        else:
            RTEC += prev
            count = 1
        prev = char
    if count > 1:
        RTEC += (str(count) + prev)
    else:
        RTEC += prev
    return RTEC

def match(read, index, distance_left, CIGAR, lower, upper, buckets, ranks, D_table, operation):
    alphabet = "actg"
    output = []
    if operation == "M":
        for char in alphabet:
            upper_new = buckets[char] + ranks[char][upper] # Will break if read char is not in genome
            lower_new = buckets[char] + ranks[char][lower]
            if upper_new == lower_new:
                continue
            if char != read[index - 1]:
                edit_distance = distance_left - 1
            else:
                edit_distance = distance_left
            if D_table[index - 2] > edit_distance:
                continue
            output.append((lower_new, upper_new, index - 1, "M" + CIGAR, edit_distance))

    if operation == "I" and distance_left > 0 and D_table[index - 2] > (distance_left - 1):
        output.append((lower, upper, index - 1, "I" + CIGAR, distance_left - 1))

    if operation == "D":
        for char in alphabet:
            upper_new = buckets[char] + ranks[char][upper] # Will break if read char is not in genome
            lower_new = buckets[char] + ranks[char][lower]
            if upper_new == lower_new:
                continue
            edit_distance = distance_left - 1
            if D_table[index - 2] > edit_distance:
                continue
            output.append((lower_new, upper_new, index, "D" + CIGAR, edit_distance))
    return output


#Stack job structure is (lower, upper, index_into_read, CIGAR_string, Edit_distance_left)

def FM_approximate_match(read, SA, buckets, ranks, D_table, distance):
    if D_table[-1] > distance:
        return []
    stack = []
    stack.extend(match(read, len(read), distance, "", len(SA), 0, buckets, ranks, D_table, "M"))
    matches = []
    while stack:
        current = stack.pop()
        if current[2] == 0:
            for i in range(current[1], current[0]):
                matches.append((SA[i], Smoke_cigar(current[3])))
            continue
        stack.extend(match(read, current[2], current[4], current[3], current[0], current[1], buckets, ranks, D_table, "M"))
        if current[4] > 0:
            stack.extend(match(read, current[2], current[4], current[3], current[0], current[1], buckets, ranks, D_table, "I"))
            stack.extend(match(read, current[2], current[4], current[3], current[0], current[1], buckets, ranks, D_table, "D"))
    return matches

def matches_to_SAM(fasta_dict, fastq_dict, name, distance):
    read_name = []
    reference_name = []
    match_index = []
    CIGARS = []
    match_string = []

    keys = [i for i in fasta_dict]

    SA_file = open(name + r"_SA.txt", "r")
    bucket_file = open(name + r"_buckets.txt", "r")
    rank_file = open(name + r"_rank.txt", "r")
    rev_SA_file = open(name + r"_rev_SA.txt", "r")
    rev_rank_file = open(name + r"_rev_rank.txt", "r")
    for i, (SA_line, bucket_line, rank_line, rev_SA_line, rev_rank_line) in enumerate(zip(SA_file, bucket_file, rank_file, rev_SA_file, rev_rank_file)):
        SA = read_SA(SA_line)
        buckets = read_bucket(bucket_line)
        ranks = read_rank(rank_line)
        rev_SA = read_SA(rev_SA_line)
        rev_ranks = read_rank(rev_rank_line)
        for read_key in fastq_dict:
            D_table = compute_D_table(fastq_dict[read_key], rev_SA, buckets, rev_ranks)
            matches = FM_approximate_match(fastq_dict[read_key], SA, buckets, ranks, D_table, distance)
            for match_cigar in matches:
                reference_name.append(keys[i])
                read_name.append(read_key)
                match_index.append(match_cigar[0] + 1)
                CIGARS.append(match_cigar[1])
                match_string.append(fastq_dict[read_key])
    SA_file.close()
    bucket_file.close()
    rank_file.close()
    rev_SA_file.close()
    rev_rank_file.close()
    output = (read_name, reference_name, match_index, CIGARS, match_string)
    return output

def print_SAM(SAM):
    for i in range(len(SAM[0])):
        sys.stdout.write(SAM[0][i] + "\t" + SAM[1][i] + "\t" + str(SAM[2][i]) + "\t" + SAM[3][i] + "\t" + SAM[4][i] + "\n")


if __name__ == '__main__':
    main()