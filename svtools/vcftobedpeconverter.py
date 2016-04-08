from svtools.bedpe import Bedpe
import re

class VcfToBedpeConverter(object):
    '''
    This is a class to take Vcf object(s) and convert them to Bedpe lines
    '''

    def __init__(self):
        '''
        Initialize a new converter
        '''
        pass

    @staticmethod
    def parse_bnd_alt_string(alt_string):
        '''
        Parse the BND alt string and return separators and region
        '''
        # NOTE The below is ugly but intended to match things like [2:222[ and capture the brackets
        # XXX Caching the compiled version of this on init may be advantageous
        r = re.compile(r'([][])(.+?)([][])')
        result = r.findall(alt_string)
        assert result
        sep1, region, sep2 = result[0]
        assert sep1 == sep2
        chrom2, breakpoint2 = region.split(':')
        breakpoint2 = int(breakpoint2)
        return sep1, chrom2, breakpoint2

    def bnd_breakpoints(self, vcf_variant):
        '''
        Return a tuple containing calculated breakpoints and orientations for a BND variant
        '''
        chrom1 = vcf_variant.chrom
        breakpoint1 = vcf_variant.pos
        orientation1 = orientation2 = '+'
        sep, chrom2, breakpoint2 = self.parse_bnd_alt_string(vcf_variant.alt)

        if vcf_variant.alt.startswith(sep):
            orientation1 = '-'
            breakpoint1 -= 1

        if sep == '[':
            orientation2 = '-'
            breakpoint2 -= 1

        # FIXME I don't think this makes ANY sense
        vcf_variant.set_info('END', breakpoint2)

        return (chrom1,
                breakpoint1,
                breakpoint1,
                chrom2,
                breakpoint2,
                breakpoint2,
                orientation1,
                orientation2)

    @staticmethod
    def simple_breakpoints(vcf_variant):
        '''
        Return a tuple containing breakpoints and orientations for simple SVs
        '''
        breakpoint1 = vcf_variant.pos
        try:
            breakpoint2 = int(vcf_variant.info['END'])
        except KeyError:
            raise ValueError('END entry in VCF required for conversion to BEDPE')

        orientation1 = orientation2 = '+'

        if 'STRANDS' in vcf_variant.info:
            strands = vcf_variant.info['STRANDS']
            orientation1, orientation2 = strands[:2]

        return (vcf_variant.chrom,
            breakpoint1,
            breakpoint1,
            vcf_variant.chrom,
            breakpoint2,
            breakpoint2,
            orientation1,
            orientation2)

    @staticmethod
    def adjust_coordinate(vcf_variant, info_tag, start, end):
        '''
        Return adjusted start and end coordinates according to the contents 
        of the tag (if it exists)
        '''
        if info_tag in vcf_variant.info:
            span = map(int, vcf_variant.info[info_tag].split(','))
            return (start + span[0], end + span[1])
        else:
            return (start, end)

    def convert(self, primary_variant, secondary_variant=None):
        '''
        Convert the passed VCF variant(s) into a BEDPE object
        '''
        vcf_variant = primary_variant
        if primary_variant is None:
            vcf_variant = secondary_variant

        try:
            sv_type = vcf_variant.info['SVTYPE']
        except KeyError:
            raise ValueError('SVTYPE field required for conversion to BEDPE')

        parser = self.simple_breakpoints
        if sv_type == 'BND':
            parser = self.bnd_breakpoints
        
        c1, s1, e1, c2, s2, e2, o1, o2 = parser(vcf_variant)

        s1, e1 = self.adjust_coordinate(vcf_variant, 'CIPOS', s1, e1)
        s2, e2 = self.adjust_coordinate(vcf_variant, 'CIEND', s2, e2)

        info_a = vcf_variant.get_info_string()
        if primary_variant is None:
            info_a = "MISSING"
            c1, s1, e1, o1, c2, s2, e2, o2 = c2, s2, e2, o2, c1, s1, e1, o1

        info_b = '.'
        if sv_type == 'BND':
            if secondary_variant is None:
                info_b = "MISSING"
            else:
                info_b = secondary_variant.get_info_string()

        # For MANTA single-ended BNDs, EVENT is not present. 
        # XXX This has probably already been calculated outside of this method. May be a candidate to memoize or otherwise cache? 
        # By adding to the variant class, perhaps?
        name = vcf_variant.var_id
        if 'EVENT' in vcf_variant.info:
            name = vcf_variant.info['EVENT']

        return Bedpe(map(str,[
            c1,
            max(s1, 0),
            max(e1, 0),
            c2,
            max(s2, 0),
            max(e2, 0),
            name,
            vcf_variant.qual,
            o1,
            o2,
            sv_type,
            vcf_variant.filter,
            info_a,
            info_b,
            vcf_variant.get_format_string(),
            vcf_variant.get_gt_string()]))

