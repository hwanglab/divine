"""
Script to convert varant annotated vcf file to excel file with all the
fields in VCF tab and the column definitions in column definitions tab

Usage:
   python vcf2xls_varant.py -i input.vcf -o output.xls -s sample1 -s sample2
   -g /gcn/data/genescores.txt -c /gcn/data/efftab_def.tsv

"""
import xlwt
import datetime
import logging
import gcn.lib.io.vcf as vcf
from gcn.lib.utils.filter import *
from collections import defaultdict
from gcn.lib.varann.vartype.varant.varant_parser import *
URLS = {'OMIMID': '"http://omim.org/entry/%s";"%s"',
        'Gene_Name':
            '"http://www.genecards.org/cgi-bin/carddisp.pl?gene=%s";"%s"',
         'OMIM_Ids': '"http://omim.org/entry/%s";"%s"',
        }
'''
VCFREM = ['ABHet', 'ABHom', 'AC', 'AF', 'AN', 'Agilent3p', 'Agilent5p',
          'AgilentCore', 'AgilentIntervalLen',
          'BaseQRankSum', 'CLNACC', 'CLNALLE',
          'CLNDBN', 'CLNDSDB', 'CLNDSDBID',
          'CLNHGVS', 'CLNORIGIN', 'CLNSIG', 'CLNSRC',
          'CLNSRCID', 'DB', 'DS', 'Dels', 'END', 'Illumina3p',
          'Illumina5p', 'IlluminaCore', 'IlluminaIntervalLen',
          'MLEAC', 'MLEAF', 'ReagentPrimary3p', 'ReagentPrimary5p',
          'ReagentPrimaryCore', 'ReagentPrimaryIntervalLen', 'OND',
          'RPA', 'RU', 'ReadPosRankSum', 'STR', 'VC', 'VQSLOD', 'culprit',
          'set', 'NEGATIVE_TRAIN_SITE' 'POSITIVE_TRAIN_SITE', 'VARANT_GENIC',
           'VARANT_INTERGENIC', 'VARANT_IMPACTCODE']
'''           
VCFREM = ['ABHet', 'ABHom', 'AC', 'AF', 'AN', 'Agilent3p', 'Agilent5p',
          'AgilentCore', 'AgilentIntervalLen',
          'BaseQRankSum', 'CLNACC', 'CLNALLE',
          'CLNDBN', 'CLNDSDB', 'CLNDSDBID',
          'CLNHGVS', 'CLNORIGIN', 'CLNSIG', 'CLNSRC',
          'CLNSRCID', 'DB', 'DS', 'Dels', 'END', 'Illumina3p',
          'Illumina5p', 'IlluminaCore', 'IlluminaIntervalLen',
          'MLEAC', 'MLEAF', 'ReagentPrimary3p', 'ReagentPrimary5p',
          'ReagentPrimaryCore', 'ReagentPrimaryIntervalLen', 'OND',
          'RPA', 'RU', 'ReadPosRankSum', 'STR', 'VC', 'VQSLOD', 'culprit',
          'set', 'NEGATIVE_TRAIN_SITE' 'POSITIVE_TRAIN_SITE', 'VARANT_GENIC',
           'VARANT_INTERGENIC', 'VARANT_IMPACTCODE']

HYPER = xlwt.easyxf('font:underline single;font:colour blue')

IN_TRANS = ['trans_id', 'region', 'splice', 'exon', 'mutation', 'codon',
            'aa', 'sift', 'pp2']
IN_VAR = ['MIM_IDS', 'MIM_PHENS', 'GAD_PHENS']

IN_EFF = ['GWASPhenotype', 'GerpRSScore', 'DB137', 'KGDB', 'KGAF',
          'ESPDB', 'ESPAF', 'EXACDB', 'EXACAF', 'IMHQ', 'IMLQ']
IN_TRANS_HDRS = ['Transcript', 'Region', 'SpliceSite', 'Exon', 'Mutation',
                 'Codon Change', 'Amino Acid Change', 'SIFT', 'Polyphen']
IN_VAR_HDRS = ['OMIMIds', 'OMIMPhenotype', 'GAD_PHENS']

VERSION = '2.0.2'


def columnnames(parser, samples=None):
    columns = [('CHROM', 'text'),
               ('POS', 'number'),
               ('ID', 'text'),
               ('REF', 'text'),
               ('ALT', 'text'),
               ('QUAL', 'number'),
               ('FILTER', 'text')]
    td = {vcf.toint: 'number',
          vcf.tofloat: 'number',
          bool: 'boolean',
          }
    INFO = parser.meta['INFO']
    infokeys = [(key, td.get(INFO[key][1], 'text')) for key in INFO]
    infokeys.sort()
    remcols = []
    for el in infokeys:
        if el[0] in VCFREM:
            remcols.append(el)
    for el in remcols:
        if el in infokeys:
            infokeys.remove(el)
    columns.extend(infokeys)
    FORMAT = parser.meta['FORMAT']
    formkeys = [(key, td.get(FORMAT[key][1], 'text')) for key in FORMAT]
    if ('AD', 'number') in formkeys:
        formkeys.remove(('AD', 'number'))
    if ('PL', 'number') in formkeys:
        formkeys.remove(('PL', 'number'))
    samplehdrs = []
    if not samples:
        samples = parser.samples
    for s in samples:
        cols = [('%s_%s' % (s, key[0]), key[1]) for key in formkeys]
        samplehdrs.extend(cols)
        #columns.extend(cols)
    return columns, infokeys, formkeys, samplehdrs


def write_cell(row, colidx, key, value, value_type, sep=','):
    url = URLS.get(key, None)
    if url and value is not None:
        term = value[0]
        pos = term.find('(')
        if pos > 0:
            term = term[:pos]
    else:
        url = None
    if isinstance(value, list):
        value = sep.join(str(el) for el in value)
        value_type = 'text'

    if value is None or value == '.':
        row.set_cell_blank(colidx)
    elif url:
        value = xlwt.Formula('HYPERLINK(' + url % (term, value) + ')')
        row.write(colidx, value, HYPER)
    else:
        try:
            getattr(row, 'set_cell_' + value_type)(colidx, value)
        except ValueError:
            row.set_cell_text(colidx, value)


def write_varant_links(row, rec, vadict, samples, formkeys, cols,
                       tmpidx, genescores=[]):
    if not isinstance(vadict, dict):
        return
    Gene_Name = vadict['GENE']
    Gene = ''
    if ':' in Gene_Name:
        gen = Gene_Name.split(':')
        genes = [i.split('(')[0] for i in gen]
        if 'NONE' not in genes:
            dist = [int(i.split('=')[-1].strip(')')) for i in gen]
            mind = min(dist)
            Gene = genes[dist.index(mind)]
        elif genes[1] == 'NONE' and genes[0] == 'NONE':
            Gene = None
        else:
            genes.remove('NONE')

            Gene = genes[0]

    else:
        Gene = Gene_Name

    if 'VARANT_IMPACTCODE' in rec.info:
        row.set_cell_text(tmpidx, rec.info['VARANT_IMPACTCODE'])
        tmpidx = tmpidx + 1
    else:
        pass
    if genescores:
        score = [i[1] for i in genescores if i[0] == Gene_Name]
        if score:
            row.set_cell_text(tmpidx, score[0])
        else:
            row.set_cell_text(tmpidx, 'A')
        tmpidx = tmpidx + 1
    else:
         pass
    url = URLS['Gene_Name']
    if Gene is not None:
        row.write(tmpidx, xlwt.Formula('HYPERLINK(' + \
                        url % (Gene, Gene_Name) + ')'), HYPER)
    else:
        row.set_cell_text(tmpidx, Gene_Name)
    tmpidx = tmpidx + 1
    for field in IN_TRANS:
        if field in vadict['TRANSCRIPT'][0]._fields:
            temp = getattr(vadict['TRANSCRIPT'][0], field)
            if field in ['sift', 'pp2']:
                eachf = []
                indiv = temp.split('__')
                for ind in indiv:
                    eachf.append(ind.split('_')[0])
                row.set_cell_text(tmpidx, ','.join(eachf))
            else:
                row.set_cell_text(tmpidx, temp)
            tmpidx += 1

    for eff in IN_VAR:
        e = vadict[eff]
        if isinstance(e, list):
            f = ','.join(e)
            write_cell(row, tmpidx, eff, f, 'text')
            tmpidx += 1

    for el in IN_EFF:
        vals = [i[0] for i in cols]
        types = [i[1] for i in cols]
        if el in vals:
            t = types[vals.index(el)]
            v = getattr(rec.info, el)
            if v:
                write_cell(row, tmpidx, el, v, t)
            else:
                write_cell(row, tmpidx, el, None, t)
            tmpidx += 1
        else:
            write_cell(row, tmpidx, el, '.', 'text')
            tmpidx += 1
    return tmpidx


def write_links(sheet, key, rec, values):
    if not isinstance(values, list):
        return
    nrow = len(sheet.rows)
    url = URLS[key]
    for v in values:
        row = sheet.row(nrow)
        row.set_cell_text(0, rec.chrom)
        row.set_cell_number(1, rec.pos)
        row.set_cell_text(2, rec.ref)
        row.set_cell_text(3, ';'.join(rec.alt))
        if rec.qual == '.':
            row.set_cell_blank(4)
        else:
            row.set_cell_number(4, rec.qual)
        row.set_cell_text(5, ';'.join(rec.filter))
        row.write(6, xlwt.Formula('HYPERLINK(' +
                                    url % (v, v) + ')'), HYPER)
        nrow += 1


def create_book(colnames, samplehdrs, opthdrs):
    logger.info('Creating book with 3 sheets(Column Def,\
                  VCF, Inheritance Rules)')
    book = xlwt.Workbook()
    sheets = {}
    col_def = book.add_sheet('Column Def')
    sheet1 = book.add_sheet('VCF')
    inh_rul = book.add_sheet('Inheritance Rules')
    col_def.row(0).set_cell_text(0, 'Column Name')
    col_def.row(0).set_cell_text(1, 'Definition')
    row = sheet1.row(0)
    idx = 0
    # col names for main sheet
    for cname, t in colnames:
        if cname not in IN_EFF:
            row.set_cell_text(idx, cname)
            idx += 1
    for hd in opthdrs:
        row.set_cell_text(idx, hd)
        idx += 1
    row.set_cell_text(idx, 'Gene_Name')
    idx = idx + 1
    for el in IN_TRANS_HDRS:
        row.set_cell_text(idx, el)
        idx += 1
    for el in ['OMIM_IDs', 'OMIM_Phenotype', 'GAD_PHENS']:
        row.set_cell_text(idx, el)
        idx += 1
    for el in IN_EFF:
        row.set_cell_text(idx, el)
        idx += 1
    for cname, t in samplehdrs:
        row.set_cell_text(idx, cname)
        idx += 1

    sheets['vcf'] = sheet1
    sheets['coldef'] = col_def
    sheets['inhrul'] = inh_rul
    return book, sheets


def write(book, outfile):
    logger.info('Saving the created workbook')
    book.save(outfile)


def run(infile, outfile, hdrfile, genefile, samples=None, fl=[[], [], [], []]):
    logger.info('Running vcf to xls conversion script')
    genescores = []
    if genefile:
        g = open(genefile, 'r')
        genescores = [(i.strip('\n').split('\t')) for i in g.readlines()]
    if samples:
        parser = vcf.VCFParser(infile, samples)
    else:
        parser = vcf.VCFParser(infile)
        samples = parser.samples

    # Write the column definitions to a tab
    col_map = defaultdict(str)
    ikeys = parser.meta['INFO'].keys()
    hf = open(hdrfile, 'r')
    col_map = {i.strip('\n').split('\t')[0]: i.strip('\n').split('\t')[1]
              for i in hf}
    for key in ikeys:
        col_map[key] = parser.meta['INFO'][key][-1]
    colnames, infokeys, formkeys, samplehdrs = columnnames(parser, samples)
    opthdrs = []
    if 'VARANT_IMPACT_CODE' in infokeys:
        opthdrs.append('VARANT_IMPACT_CODE')
    if genefile:
        opthdrs.append('Project_Gene_score')
    book, sheets = create_book(colnames, samplehdrs, opthdrs)
    tempidx = 1
    logger.info('Writing the column definitions in Column Def sheet')
    for key, value in sorted(col_map.items()):
        if [j[0] for j in IN_TRANS_HDRS + opthdrs +  ['Gene_Name'] \
         + IN_VAR_HDRS if key.lower() in j.lower()] or\
                    [j[0] for j in colnames if key.lower() in j[0].lower()]:
            sheets['coldef'].row(tempidx).set_cell_text(0, key)
            sheets['coldef'].row(tempidx).set_cell_text(1, value)
            tempidx = tempidx + 1
    f = Filter(fl)
    rowidx = 1
    senum = 1
    excelrows = 0
    basename = outfile[:-4]
    booknum = 2
    logger.info('Writing data to VCF tab')
    logger.info('The fields from the varant annotation are '\
                 + ','.join(IN_TRANS_HDRS + IN_VAR_HDRS))
    for rec in parser:
        parser.parseinfo(rec)
        if excelrows > 60000:
            logger.info('Write out book and create a new one as the number\
                         of rows reached 60000')
            write(book, outfile)
            outfile = basename + '_' + str(booknum) + '.xls'
            book, sheets = create_book(colnames, samplehdrs)
            tempidx = 1
            for key, value in sorted(col_map.items()):
                if [j[0] for j in IN_TRANS_HDRS + opthdrs +  ['Gene_Name'] \
                               + IN_VAR_HDRS if key.lower() in j.lower()] or\
                               [j[0] for j in colnames if key.lower() in j[0].lower()]:
                    sheets['coldef'].row(tempidx).set_cell_text(0, key)
                    sheets['coldef'].row(tempidx).set_cell_text(1, value)
                    tempidx = tempidx + 1
            excelrows = 0
            booknum += 1
            rowidx = 1
            senum = 1
        if f.retain(rec):
            pass
        else:
            continue
        excelrows += 1
        parser.parsegenotypes(rec)
        row = sheets['vcf'].row(rowidx)
        idx = 0
        for e, t in colnames[:7]:
            if e is 'FILTER' and rec[e.lower()][0] is '.':
                write_cell(row, idx, e.lower(), None, t, ';')
            else:
                write_cell(row, idx, e.lower(), rec[e.lower()], t, ';')
            idx += 1
        info = rec.info
        for e, t in infokeys:
            if e not in IN_EFF:
                write_cell(row, idx, e, info.get(e, None), t)
                idx += 1
        par_ant = parse(info)
        ga = get_prior_geneannot(info, alltrans=False)
        idx = write_varant_links(row, rec, ga, samples, formkeys,
                                 colnames, idx, genescores)
        senum += 1
        sheets['vcf'].panes_frozen = True
        sheets['vcf'].remove_splits = True
        sheets['vcf'].vert_split_pos = 2
        sheets['vcf'].horz_split_pos = 1
        for s in samples:
            ss = getattr(rec, s)
            for e, t in formkeys:
                v = ss.get(e, None)
                write_cell(row, idx, e, v, t)
                idx += 1

        rowidx += 1
        if not excelrows % 500:
            sheets['vcf'].flush_row_data()
    write(book, outfile)


def get_logger():
    #create log file with timestamp
    now = str(datetime.datetime.now())
    logfile = "vcf2xls_varant_" + now.replace(' ', '.') + ".log"
    logger = logging.getLogger('vcf2xls_varant')
    logger.setLevel(logging.DEBUG)
    # create file handler which logs even debug messages
    fh = logging.FileHandler(logfile)

    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - \
                                %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    fh.setLevel(logging.DEBUG)
    logger.addHandler(fh)
    return logger

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Filter variants \
                                                     in the vcf')
    parser.add_argument('-i', dest='infile', help='input vcf file')
    parser.add_argument('-s', dest='samples', action='append',
                          help='Samples to be displayed in output\
                                (order is as in the vcf), \
                               format is -s sample1 -s sample2')
    parser.add_argument('-o', dest='outfile', default=None,
                         help='output file')
    parser.add_argument('-f', dest='filterlist', default="",
                         help='filterfile')
    parser.add_argument('-g', dest='genefile', default=None,
                         help='File with genelist and the corresponding score')
    parser.add_argument('-c', dest='hdrfile', default=None,
                         help='Column Definitions file')
    options = parser.parse_args()

    fl = []
    logger = get_logger()
    logger.info('software version: ' + VERSION)
    logger.info('Input VCF' + options.infile)
    if options.samples:
        logger.info('Samples for which the excel is to be generated are  ' + \
                    ','.join(options.samples))
    if options.filterlist:
        logger.info('Filter file ' + options.filterlist)
    logger.info('File with project gene score ' + options.genefile)
    logger.info('Column Definitions file ' + options.hdrfile)
    if options.filterlist:
        for el in open(options.filterlist, 'rU'):
            fl.append(el.strip().split(','))
    run(options.infile, options.outfile, options.hdrfile,
         options.genefile, options.samples, [[], [], [], []])
