import unittest
import NameSplitter
import re

class TestNameSplitter(unittest.TestCase):

    def setUp(self):
        self.splitter = NameSplitter.NameSplitter()

    def test_split(self):
        tests = '''
            # All test cases are taken from real crowd-sourced data.
            
            ### Basics: split on punct or preposition, then reconnect 'jr' and 'sr'.
            
            Stuttgart | Stuttgart
            James D. Ray Jr. | James D. Ray Jr.
            C. Earle Smith, Jr. | C. Earle Smith, Jr.
            Olga Lakela, Jackie Patman | Olga Lakela | Jackie Patman
            James D. Ray Jr., C. Earle Smith, Jr., Olga Lakela, Jackie Patman | James D. Ray Jr. | C. Earle Smith, Jr. | Olga Lakela | Jackie Patman
            D. S. Correll and Helen B. Correll | D. S. Correll | Helen B. Correll
            R. K Godfrey with Angus Gholson | R. K Godfrey | Angus Gholson
            H. Maurushat; V. Sullivan, C. Hudson; D. Wise, R.K. Godfrey | H. Maurushat | V. Sullivan | C. Hudson | D. Wise | R.K. Godfrey
            Collected by: R. K. Godfrey | R. K. Godfrey
            william p adams | william p adams
            A. Gholson Jr w/Wilson Baker | A. Gholson Jr | Wilson Baker
            D. B. Ward, with H. F. Decker | D. B. Ward | H. F. Decker
            Cecil R. Slaughter, Ph.D. | Cecil R. Slaughter, Ph.D.
            George Eiten, Liene T. Eiten, Gil M. Felippe & J.M. de Freitas Campes | George Eiten | Liene T. Eiten | Gil M. Felippe | J.M. de Freitas Campes
            Tim Reeves, Sr | Tim Reeves, Sr
            John W/ Thieret | John | Thieret
            A. Gholson, Jr., w/ Dr. Bob Godfrey & D. C. Vickers | A. Gholson, Jr. | Dr. Bob Godfrey | D. C. Vickers
            Gholson, Jr., with Dr Bob Godfrey and D.C. Vickens | Gholson, Jr. | Dr Bob Godfrey | D.C. Vickens
            Charles T. Bryson (and) Will McDearman | Charles T. Bryson | Will McDearman
            S Leonard + H McAninch | S Leonard | H McAninch
            Mrs. C. W. Hanes | Mrs. C. W. Hanes
            Coll. Ame Garthwright | Coll. Ame Garthwright
            W.H. Wagner , Jr. | W.H. Wagner , Jr.
            Ann F-Johnson | Ann F-Johnson
            R.R. Snelling Collectors | R.R. Snelling
            J.R.Powers.Collr. | J.R.Powers.
            J.R Powers collr | J.R Powers
            K.S.H. | K.S.H.
            R.wittle | R.wittle
            H. Ruckes, Jr. | H. Ruckes, Jr.
            M.rhen.collr | M.rhen.
            E. P. VanDuzee | E. P. VanDuzee
            Drake & Hottes | Drake | Hottes
            Angusto-Cylin Dricus | Angusto-Cylin Dricus
            robert denoble collection | robert denoble collection # TODO: more adhoc cleaning?
            P.D .Hurd | P.D .Hurd
            e.e.gilbert & c.d.mac neil | e.e.gilbert | c.d.mac neil
            
            
            ### Commas which should be periods:
            
            # 98/12548 in herbarium match m{\b[A-Z],}
            
            R, Kral & P.L. Redfearn | R. Kral | P.L. Redfearn
            R, K, Godfrey | R. K. Godfrey
            W. A, Sliveus | W. A. Sliveus
            R,K, Godfrey & J.P, Gillespie | R.K. Godfrey | J.P. Gillespie
            R. K,. Godfrey and Richard D. Houk | R. K. Godfrey | Richard D. Houk
            J. D, McCarty | J. D. McCarty
            R,O, Schuster | R.O. Schuster
            
            
            ### Name distribution:
            
            # TODO: frequency?
            
            Nancy Craft Coile, w/ Robert, Danielle & Robbie Coile | Nancy Craft Coile | Robert Coile | Danielle Coile | Robbie Coile
            P.J. Crutchfeld & Laura & Thomas Crutchfield | P.J. Crutchfeld | Laura Crutchfield | Thomas Crutchfield
            D. B. & S. S. Ward | D. B. Ward | S. S. Ward
            Robert & Mabel Kral/ | Robert Kral | Mabel Kral
            R. K. Godfrrey with Robt. & John Lazor | R. K. Godfrrey | Robt. Lazor | John Lazor
            Bruce Hansen with T.&B. Cochrane, C.S. Keller & M. Waterway | Bruce Hansen | T. Cochrane | B. Cochrane | C.S. Keller | M. Waterway
            R & JK Robertson | R Robertson | JK Robertson
            J. A. Chemsak, A. & M. Michchelbacher & W.W. Middlekauff | J. A. Chemsak | A. Michchelbacher | M. Michchelbacher | W.W. Middlekauff
            J.A. & M.A. CHEMSAK, E.G. & J.M. LINSLEY | J.A. CHEMSAK | M.A. CHEMSAK | E.G. LINSLEY | J.M. LINSLEY
            D. Spencer; R., J. &A. Ryckman | D. Spencer | R. Ryckman | J. Ryckman | A. Ryckman
            JM and SM Burns | JM Burns | SM Burns
            r & j. robertson | r robertson | j. robertson
            J&R. Robertson | J Robertson | R. Robertson
            J.A & M.A Chemsak + E.G & J.M Linsley | J.A Chemsak | M.A Chemsak | E.G Linsley | J.M Linsley
            
            # TODO: Is this correct? "Lafon" and "Grey" are weird first names.
            # If a name list look-up were incorporated, would the behavior be different?
            
            Lafon & Gray Bill | Lafon Bill | Gray Bill
            
            
            ### Slashes:
            
            # 242/12548 in herbarium match m{/}
            
            Lytton J. Musselman / Elizabeth R. Musselman | Lytton J. Musselman | Elizabeth R. Musselman
            Michel G. Lelong / Ken Rogers | Michel G. Lelong | Ken Rogers
            SW Leonard/D. Culwell/M. Ripperton | SW Leonard | D. Culwell | M. Ripperton
            A. H. Curtiss/ det. C. B. Heiser, Jr. | A. H. Curtiss | det. C. B. Heiser, Jr.
            # R/K/ Godfrey & John Morrill
            
            
            ### ACK!! No punctuation between names:
            
            # 94/12548 in herbarium match m{^[A-Za-z. ]+$} && /(.*\.){4,}/ && ! /\band|with\b/
            
            # M.B. H.L. | M.B. | H.L.
            # D.R.Windler B.r. Sinor | D.R.Windler | B.r. Sinor
            # S.W.Leonard D. Culwell M.Ripperton | S.W.Leonard | D. Culwell | M.Ripperton
            # R. K. Godfrey Richard D. Houk | R. K. Godfrey | Richard D. Houk
            # R.&A.Ryckman C.Christianson | R. Ryckman | A.Ryckman | C.Christianson
            
            
            ### Parens:
            
            # 287/12548 in herbarium match m{[()]}
            # 175/12548 in herbarium match m{[()]} && ! m{\(\?\)}
            
            # TODO: perhaps an earlier phase in the process should remove parenthetical expressions
            # which include dates? Are these determinations rather than the original collection?
            
            John Mayberg (By W) | John Mayberg | By W
            (Mary L. Leigh) J. Rowntrey | Mary L. Leigh | J. Rowntrey
            A.R. Diamond (w. J.D. Freeman) | A.R. Diamond | w. J.D. Freeman # TODO: special handling for "w."?
            David Hall (w/ Gary Schultz) | David Hall | Gary Schultz
            R K Godfrey ( Shirley Mah Kooyman 1980) | R K Godfrey | Shirley Mah Kooyman 1980
            Loran C Anderson ( Scott Sundberd 1987) | Loran C Anderson | Scott Sundberd 1987
            (MARY L. LEIGH) J. ROUNTREY | MARY L. LEIGH | J. ROUNTREY
            
            # (Karl, Godfrey 1958); R K Godfrey 1976 | ???
            # R. K. Godfrey (det.) & Richard D. Houk | R. K. Godfrey (det.) | Richard D. Houk # TODO: maybe strip out the "(det.)"?
            
            ### Dashes:
            
            # 226/12548 in herbarium match m{-}
            # 128/12548 in herbarium match m{-} && m{^[A-Z -]+$} (Perhaps all from a single user?)
            
            Barton H. warnock, Reginald Rose-Innes | Barton H. warnock | Reginald Rose-Innes
            Marie-Victorin, Rolland-Germain, Marcel Raymond | Marie-Victorin | Rolland-Germain | Marcel Raymond
            REGINALD TOSE-INNES & BARTON H. WARNOCK | REGINALD TOSE-INNES | BARTON H. WARNOCK
            
            # TODO: I really have no good heuristic in mind for these:
            
            # H.H. Iltis - D. Parker
            # ROBERT K GODFREY - R W SIMONS - ANGUS GHOLSON
            # Robert F- Thorne
            # Steve L. Orzell- Edwin L. Bridges
            # A GHOLSON, JR - SUSANNE COOPER - WILSON BAKER
            # H Maurushat - V. Sullivan, C. Hudson
            # R.K. Godfrey w/- Christopher Campbell
            
            
            ### Inverted Names:
            
            # rare?
            
            # Betts, Thealcald Jonas, Baker.
            Powers, J. R. | J. R. Powers
            HAGEN, K. S. |  K. S. HAGEN
            Stage, G., Snelling, R.R. | G. Stage | R.R. Snelling
            Powell, J. | J. Powell
            
            
            ### Question marks:
            
            # 323/12548 in herbarium match m{\?}
            # TODO: Just drop, perhaps with whole phrase?
            # Maybe this is resolved at an earlier step in the processing?
            
            # Judith Canne and Jose Schunkell (or Schunkey?)
            # ?Rufus Crane
            # Lloyd T. (?Y.?) Card (?Cart?)
            
            
            ### Just initials:
            
            # 122/12548 in herbarium match !m{\w{2}}
            
            # A.H.S.F.
            
            
            ### Dates:
            
            # 77/12548 in herbarium match m{\b(19|20)\d\d\b}
            # TODO: Dates removed at an earlier step in processing?
            # These records often also involve a determination.
            
            # Donald Eves (1956) and L. J. Ultal (1981)
            
            
            ### Other numbers:
            
            # 47/12548 in herbarium match m{\b\d{3}\b}
            
            # D.B. Ward 3-19, with BTY 421
            
            
            ### Other data that shouldn't have been included:
            
            # J.R. Powers Collr, UC Berkely EMEC
            # JS Buckett M.R & RC Gardener Coll. det. W.D.Sumlin 19173
            
            
            ### Random punctuation:
            
            # 30/12548 in herbarium match m{[^A-Za-z0-9();/+&,. -?]} && !m{\?}
            
            # R> K> Godfrey with Robt. & John Lazor | R> K> Godfrey | Robt. Lazor | John Lazor
            # Loran C. Anderson w/Gil Nelson R>K> Godfrey Herbarium (FSU)
            # R:D. Houk ans R:K: Godfrey
            # Loran C: Anderson
        '''
        no_comments = [re.match(r'^[^#]*', line).group() for line in tests.split('\n')]
        arrays = [
            [element.strip(' ') for element in line.split('|')]
                for line in no_comments]
        fail = 0;
        for (input, expected) in [(an_array[0],an_array[1:]) for an_array in arrays if an_array != ['']]:
            actual = self.splitter.split(input)
            if actual != expected:
                print('FAIL: given "' + input + '" expected ', expected, '; actual ', actual)
                fail += 1
        self.assertFalse(fail)        

if __name__ == '__main__':
    unittest.main()