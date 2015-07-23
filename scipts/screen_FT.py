import sys, os, glob
import argparse
from operator import itemgetter
from subprocess import Popen, PIPE

def main(args):
    
    vcfs = reduce(lambda x,y: x.extend(y),args.vcf)
    #check if vcftools are available
    if which("vcftools") is not None:
    #build the command
        cl = ["vcftools",  "--recode", "--recode-INFO-all", "--stdout"]
        if args.PASS:
            cl.append("--keep-filtered")
            cl.append("PASS")
    else:
        print "vcftools not available on commadn line, please install vcftools and add it to your $PATH variable"
        return 1
    
    for vcf in vcfs:
        #vcf_base_name = ".".join(os.path.split(vcf)[-1].split(".")[0:-1]) # collect the base name without the final .vcf
        cl.append("--vcf")
        cl.append("{}".format(vcf))
        p = Popen(cl, stdin=PIPE, stdout=PIPE, stderr=None)
        for line in p.stdout:
            if line.startswith("#") or line.startswith("="):
                continue
            variation   = line.rstrip().split("\t")[7]
            description = dict(item.split("=") for item in variation.split(";"))
            passed = True
            if args.max_cov:
                if float(description["COVA"]) > args.max_cov or float(description["COVB"]) > args.max_cov:
                    passed = False and passed
            if args.min_cov:
                if float(description["COVA"]) < args.min_cov or float(description["COVB"]) < args.min_cov:
                    passed = False and passed
            if args.min_number_links:
                if float(description["LTE"]) < args.min_number_links:
                    passed = False and passed
            if args.min_ratio_links:
                if float(description["LTE"])/float(description["LFW"]) < args.min_ratio_links:
                    passed = False and passed

            if passed:
                print line

    
    return 0


# vcftools  P1426/P1426_104_intra_chr_events.vcf


def which(program):
    #taken from http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


if __name__ == '__main__':
    parser = argparse.ArgumentParser("""
    This script takes as input a vcf file produced by FindTransloations and applies filters specifed as option. vcftools needs to be available on $PATH"
    """)
    parser.add_argument('--vcf', type=str, required=True, action='append', nargs='+', help="vcf file geenrated by find translocations")
    parser.add_argument('--keep-filtered', default=False,  help="keep only PASS entries", action='store_true', dest='PASS')
    parser.add_argument('--max-cov',          type=float,  help="max coverage on either window A and window B")
    parser.add_argument('--min-cov',          type=float,  help="min coverage on either window A and window B")
    parser.add_argument('--min-number-links', type=int,    help="minimum number of proper link between  window A and window B")
    parser.add_argument('--min-ratio-links',  type=float,    help="min ration between the links leaving A and the number of link falling properly on window B")
    args = parser.parse_args()

    main(args)



