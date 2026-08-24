"""Check the DerSimonian-Laird estimator against metafor, the reference implementation.

illumeta_meta.py pools cohorts in pure Python with no numpy or scipy, and every candidate
candidate CpG this tool reports comes out of that code. Until 2026-08-08 nothing had compared it
with an established implementation.

The comparison itself was run separately: 7,713 CpGs drawn from the four re-run cohorts and
stratified to force the boundary cases, scored against metafor 5.0.1
rma(yi, sei, method="DL", test="z"). Worst relative disagreement across tau2, Q, I2 and the
fixed and random effect, SE and p was 4.8e-12, with no case exceeding a 1e-9 relative
tolerance and no NA mismatches; p-values agreed to |dlog10| < 1.5e-14 as deep as 1e-25.

Those metafor outputs are frozen below so the check runs in CI, which has no R and no
metafor. The cases are chosen to exercise different code paths rather than to be
representative -- a uniform sample of an array is overwhelmingly null and homogeneous, and
would leave the interesting branches untested.

test="z" matches how illumeta_meta.py refers its Wald statistic to a standard normal. That
the normal reference is anti-conservative for k <= 5 is a stated limitation of the
analysis; it is not what this test is about, which is whether the arithmetic is right.
"""

import math
import unittest

from illumeta_meta import _random_effect_meta_one


#: Inputs and metafor 5.0.1 outputs. Values are full double precision.
METAFOR_DL_CASES = (
    # tau2 truncated to zero (Q <= df): the random model must collapse onto the fixed one
    dict(
        cpg='cg12377401',
        yi=[0.107242542869259, 0.0492312567232336],
        sei=[0.06971310941549395, 0.03699845214632537],
        tau2=0.0,
        Q=0.540281860586634,
        I2=0.0,
        fixed_effect=0.0619802239165895,
        fixed_se=0.0326810279752114,
        fixed_p=0.0578913149183489,
        random_effect=0.0619802239165895,
        random_se=0.0326810279752114,
        random_p=0.0578913149183489,
    ),
    dict(
        cpg='cg23271057',
        yi=[0.0462788566896002, 0.0472164722806028, -0.0225923272278338, 0.0462448214753404],
        sei=[0.062760811783919, 0.045842558959043, 0.0388667930998833, 0.04608463375266519],
        tau2=0.0,
        Q=2.04532099007128,
        I2=0.0,
        fixed_effect=0.0220320078093895,
        fixed_se=0.0231710649615323,
        fixed_p=0.341684893746014,
        random_effect=0.0220320078093895,
        random_se=0.0231710649615323,
        random_p=0.341684893746014,
    ),
    dict(
        cpg='cg00234998',
        yi=[0.0184720158519924, -0.00677479419597124, -0.046006361290404, -0.00616421541413503],
        sei=[0.04844974151567773, 0.07364650856826342, 0.046771050288699254, 0.044059153082323145],
        tau2=0.0,
        Q=0.94554359952048,
        I2=0.0,
        fixed_effect=-0.0111119978843877,
        fixed_se=0.0251365888777758,
        fixed_p=0.658442408868225,
        random_effect=-0.0111119978843877,
        random_se=0.0251365888777758,
        random_p=0.658442408868225,
    ),
    dict(
        cpg='cg26376241',
        yi=[-0.0024614560747136, 0.0186995425930752, 0.0408853201789148],
        sei=[0.07744466432283098, 0.04944430837136968, 0.05296860592613592],
        tau2=0.0,
        Q=0.229520611188083,
        I2=0.0,
        fixed_effect=0.0233973556691533,
        fixed_se=0.0327526900185368,
        fixed_p=0.475001881078834,
        random_effect=0.0233973556691533,
        random_se=0.0327526900185368,
        random_p=0.475001881078834,
    ),
    # high heterogeneity: tau2 dominates the weights
    dict(
        cpg='cg23712156',
        yi=[-0.198560935564223, -0.10491671232307, -0.00927596108629825, 0.0647861370214802],
        sei=[0.058816207712280744, 0.05875613575251253, 0.030858687095214395, 0.049594449970245334],
        tau2=0.00818885133194459,
        Q=13.8924950024244,
        I2=78.4056067720271,
        fixed_effect=-0.0349752218731058,
        fixed_se=0.0221651650289259,
        fixed_p=0.114580209572017,
        random_effect=-0.0559639027965151,
        random_se=0.051635739368092,
        random_p=0.278444136879069,
    ),
    dict(
        cpg='cg23119809',
        yi=[0.0691803432267499, -0.013383751633735, -0.00317593935701721, 0.153097460466095],
        sei=[0.042320670355790464, 0.04785389103937617, 0.024622775431934626, 0.03197967996692315],
        tau2=0.005676648556467,
        Q=17.0127805570796,
        I2=82.3661982241251,
        fixed_effect=0.0489319285438169,
        fixed_se=0.0166154655854242,
        fixed_p=0.00322993463507405,
        random_effect=0.0527082595451962,
        random_se=0.0419641177860984,
        random_p=0.209104500045962,
    ),
    dict(
        cpg='cg07610169',
        yi=[0.254081835920043, -0.0135592436726286, 0.130090975119004, 0.306579272840777],
        sei=[0.07009589529583973, 0.05818802250104813, 0.046740460698855395, 0.05261320145240294],
        tau2=0.016452021942169,
        Q=18.8305147447125,
        I2=84.0684121455449,
        fixed_effect=0.165424963684127,
        fixed_se=0.0275464542610245,
        fixed_p=1.90968111529009e-09,
        random_effect=0.168260488933454,
        random_se=0.0702036330484687,
        random_p=0.0165412493113143,
    ),
    dict(
        cpg='cg24607948',
        yi=[0.0557513552566489, -0.084420194312377],
        sei=[0.03994925180398519, 0.045742647872629216],
        tau2=0.00797986537724216,
        Q=5.32708561942274,
        I2=81.2280096202328,
        fixed_effect=-0.00490092074134309,
        fixed_se=0.0300894960984485,
        fixed_p=0.870614393627385,
        random_effect=-0.0125635640431829,
        random_se=0.0700633991324846,
        random_p=0.857688737478057,
    ),
    # k = 2, the smallest analysable set
    dict(
        cpg='cg17593404',
        yi=[0.128018668759075, 0.0319098578960744],
        sei=[0.05770185414297931, 0.04744562415066968],
        tau2=0.00182815615145788,
        Q=1.65518368163554,
        I2=39.5837446263446,
        fixed_effect=0.0706780296215719,
        fixed_se=0.0366476243371067,
        fixed_p=0.0537825008606087,
        random_effect=0.0743538686571615,
        random_se=0.0477257724192831,
        random_p=0.119247577379379,
    ),
    dict(
        cpg='cg09306214',
        yi=[0.175133526703474, 0.0210682823723234],
        sei=[0.06754934123901948, 0.05759466008107738],
        tau2=0.0079280205695684,
        Q=3.01217305649928,
        I2=66.801376240906,
        fixed_effect=0.0859227175491827,
        fixed_se=0.0438266847807641,
        fixed_p=0.0499360404210145,
        random_effect=0.0940579140588359,
        random_se=0.076926452574459,
        random_p=0.221443319256648,
    ),
    dict(
        cpg='cg12924880',
        yi=[-0.0214491092744944, 0.0655815441124501],
        sei=[0.028990228562797697, 0.04480190254629056],
        tau2=0.00236334540253395,
        Q=2.65986025549014,
        I2=62.4040399139041,
        fixed_effect=0.0042364989992029,
        fixed_se=0.0243391418882372,
        fixed_p=0.861817405775922,
        random_effect=0.015362963598416,
        random_se=0.0429959305685636,
        random_p=0.720858218642995,
    ),
    dict(
        cpg='cg21196708',
        yi=[-0.0684239345244859, -0.00307763234470571],
        sei=[0.04552411165071437, 0.025434949163972164],
        tau2=0.000775378914005295,
        Q=1.57026125099486,
        I2=36.3163295683163,
        fixed_effect=-0.0186234169131218,
        fixed_se=0.0222043106108747,
        fixed_p=0.401620953961418,
        random_effect=-0.0248434477854338,
        random_se=0.0307987797027627,
        random_p=0.419875415930964,
    ),
    # k = 4, every cohort contributing
    dict(
        cpg='cg05691871',
        yi=[0.118555305775358, -0.0341885957379178, -0.0439578504641851, -0.0906454571119031],
        sei=[0.06524926756208439, 0.04160361612100095, 0.04060050284319525, 0.05276217624001491],
        tau2=0.00283543125952769,
        Q=6.6302219565809,
        I2=54.7526459951719,
        fixed_effect=-0.0287513797822165,
        fixed_se=0.0237123176251372,
        fixed_p=0.225317902053407,
        random_effect=-0.0219554388993887,
        random_se=0.0362508181363778,
        random_p=0.54474479443998,
    ),
    dict(
        cpg='cg14273031',
        yi=[0.0230103652599278, -0.0189609999362808, 0.0711360339268747, -0.114220931198315],
        sei=[0.06698337317047591, 0.04200771198112256, 0.05371816697545095, 0.06757094081092734],
        tau2=0.00199104694986325,
        Q=4.90637213863961,
        I2=38.855025358272,
        fixed_effect=-0.00441373723432438,
        fixed_se=0.0271650463406707,
        fixed_p=0.870929021514367,
        random_effect=-0.00556342823387836,
        random_se=0.0358470369535184,
        random_p=0.876664349848511,
    ),
    dict(
        cpg='cg04214629',
        yi=[0.0618796194491953, -0.000854068184580559, 0.0304374241177783, 0.175619232907809],
        sei=[0.09963520881482564, 0.05312569604317248, 0.07179619295428412, 0.08051497521636855],
        tau2=0.000783391305618116,
        Q=3.43707111897897,
        I2=12.7163827529063,
        fixed_effect=0.048457490467665,
        fixed_se=0.0352825300114299,
        fixed_p=0.169623945920313,
        random_effect=0.0514770342951517,
        random_se=0.038440708506818,
        random_p=0.180528962625614,
    ),
    dict(
        cpg='cg17254807',
        yi=[-0.102236678670537, -0.117247548307053, -0.0889534228583493, 0.0387334779361748],
        sei=[0.05172856004932325, 0.05196040789304669, 0.039372139251298906, 0.05518359306978708],
        tau2=0.00184441063741589,
        Q=5.30794883083702,
        I2=43.4809924584951,
        fixed_effect=-0.0735326973776411,
        fixed_se=0.0241291555818938,
        fixed_p=0.00230782300963731,
        random_effect=-0.0709300878595928,
        random_se=0.0326097386505671,
        random_p=0.0296211408781117,
    ),
    # very small SE: inverse-variance weights are most lopsided here
    dict(
        cpg='cg09274827',
        yi=[-0.105132242297844, -0.0733033856216088, -0.0515787448090225, 0.00791531361272746],
        sei=[0.035297192335700917, 0.029428217522754662, 0.01991010952608457, 0.027391091024077746],
        tau2=0.00111558338991652,
        Q=7.54291961424394,
        I2=60.2276021298856,
        fixed_effect=-0.0496473267204729,
        fixed_se=0.0131161546435386,
        fixed_p=0.000153582271894761,
        random_effect=-0.0524046810873772,
        random_se=0.0216908477215713,
        random_p=0.0156928700249517,
    ),
    dict(
        cpg='cg18548043',
        yi=[-0.051010880639613, 0.032294781378871, -0.0258746717241021, 0.00833044237845804],
        sei=[0.02442598807860661, 0.024979073824378223, 0.044087864498955215, 0.019632804127029998],
        tau2=0.000776543359323113,
        Q=6.45603853568407,
        I2=53.531875879949,
        fixed_effect=-0.00398264014312859,
        fixed_se=0.0125120993259776,
        fixed_p=0.750255024552367,
        random_effect=-0.00615296691069627,
        random_se=0.0193511078579346,
        random_p=0.750511803617177,
    ),
    dict(
        cpg='cg21625881',
        yi=[-0.0688926991273211, 0.010219907071098, -0.0162651416091766, 0.039765864116255],
        sei=[0.02712955739342392, 0.023850980286055295, 0.017236329651863586, 0.02914510467374756],
        tau2=0.00102187754249873,
        Q=8.47244577687164,
        I2=64.591098261266,
        fixed_effect=-0.0109103720103243,
        fixed_se=0.011425947734399,
        fixed_p=0.339640021026147,
        random_effect=-0.00957708909806571,
        random_se=0.0200399149470147,
        random_p=0.632720893205182,
    ),
    dict(
        cpg='cg06331970',
        yi=[-0.0745620446177546, -0.0595344024146338, -0.00331520753718806, -0.0187851824786811],
        sei=[0.029879430121107803, 0.022092114650201274, 0.020280586473165083, 0.01993166892986364],
        tau2=0.00050908021838834,
        Q=6.02241200496258,
        I2=50.1860716681631,
        fixed_effect=-0.0321335482082298,
        fixed_se=0.0110991336221979,
        fixed_p=0.00378989192541864,
        random_effect=-0.0351305515630657,
        random_se=0.0160004432340393,
        random_p=0.0281206829670939,
    ),
    # typical probe
    dict(
        cpg='cg27569040',
        yi=[0.102073208242938, -0.171480481436042, -0.0855444262115614],
        sei=[0.07229622891259534, 0.07621922418576596, 0.05075533674762571],
        tau2=0.011690858634065,
        Q=7.41340205254226,
        I2=73.0218328127213,
        fixed_effect=-0.057468533802379,
        fixed_se=0.0364749206991514,
        fixed_p=0.115126559394596,
        random_effect=-0.0523825785149059,
        random_se=0.073255819214216,
        random_p=0.474569671321414,
    ),
    dict(
        cpg='cg24805759',
        yi=[0.0736783722929256, -0.0127023505958004, 0.166031064493541],
        sei=[0.09696192536905979, 0.07029364957193754, 0.08467327986462345],
        tau2=0.00226397237709379,
        Q=2.65248070099732,
        I2=24.598885893948,
        fixed_effect=0.0634144574869223,
        fixed_se=0.0472337891280045,
        fixed_p=0.179412695745196,
        random_effect=0.0670254833387036,
        random_se=0.0550051812484228,
        random_p=0.223022494646238,
    ),
    dict(
        cpg='cg06896044',
        yi=[-0.0380378007506419, -0.113616336893994, 0.0719915024808015],
        sei=[0.07278731790113012, 0.05243499497020311, 0.05437981588618428],
        tau2=0.00708644407102873,
        Q=6.06544548081344,
        I2=67.0263296187145,
        fixed_effect=-0.0271253540470144,
        fixed_se=0.0335083237222853,
        fixed_p=0.418221322754778,
        random_effect=-0.0264370599536483,
        random_se=0.0595511429313336,
        random_p=0.657086867679478,
    ),
    dict(
        cpg='cg14989202',
        yi=[-0.02251419355742, 0.0218239797238251, 0.184467423968857],
        sei=[0.06586707077528375, 0.028503613328174467, 0.0566544572562377],
        tau2=0.00692326062996477,
        Q=7.77159348545526,
        I2=74.2652519879861,
        fixed_effect=0.0446411821527151,
        fixed_se=0.0237497882268366,
        fixed_p=0.060156410189155,
        random_effect=0.0600675604968831,
        random_se=0.0561873685349257,
        random_p=0.285043497914189,
    ),
)


class MetaforEquivalenceTests(unittest.TestCase):
    """_random_effect_meta_one must reproduce metafor to near machine precision."""

    #: Both implementations do the same arithmetic in IEEE doubles but accumulate it in a
    #: different order, so agreement should be near machine precision. A materially looser
    #: agreement would mean the estimators differ, not that rounding accumulated.
    REL_TOL = 1e-9

    def _assert_close(self, field, ours, reference, cpg):
        if math.isnan(reference):
            self.assertTrue(math.isnan(ours), f'{cpg}/{field}: metafor NA, ours {ours!r}')
            return
        self.assertFalse(math.isnan(ours), f'{cpg}/{field}: ours NA, metafor {reference!r}')
        scale = max(abs(ours), abs(reference))
        rel = abs(ours - reference) / scale if scale > 0 else 0.0
        self.assertLessEqual(
            rel, self.REL_TOL,
            f'{cpg}/{field}: ours={ours!r} metafor={reference!r} rel={rel:.3e}',
        )

    def test_matches_metafor_dersimonian_laird(self):
        for case in METAFOR_DL_CASES:
            with self.subTest(cpg=case['cpg']):
                yi, sei = case['yi'], case['sei']
                got = _random_effect_meta_one(yi, sei, [True] * len(yi))
                self.assertEqual(got['k'], len(yi))
                for field in ('tau2', 'Q', 'I2', 'fixed_effect', 'fixed_se',
                              'fixed_p', 'random_effect', 'random_se', 'random_p'):
                    self._assert_close(field, got[field], case[field], case['cpg'])

    def test_zero_tau2_collapses_random_onto_fixed(self):
        """When Q <= df the DL estimate truncates to zero and the two models coincide.

        Worth pinning separately: it is the branch a refactor is most likely to break,
        and breaking it would silently change every homogeneous CpG in the analysis.
        """
        checked = 0
        for case in METAFOR_DL_CASES:
            if case['tau2'] != 0.0:
                continue
            checked += 1
            yi, sei = case['yi'], case['sei']
            got = _random_effect_meta_one(yi, sei, [True] * len(yi))
            self.assertEqual(got['tau2'], 0.0, case['cpg'])
            self.assertAlmostEqual(got['random_effect'], got['fixed_effect'], places=15)
            self.assertAlmostEqual(got['random_se'], got['fixed_se'], places=15)
        self.assertGreater(checked, 0, 'fixture lost its tau2 == 0 cases')


if __name__ == '__main__':
    unittest.main()
