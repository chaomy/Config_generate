#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2018-07-23 15:14:07


from optparse import OptionParser
import gn_dd_data_hcp
import gn_dd_data_bcc

if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option('-t', "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")
    (options, args) = parser.parse_args()

    bcc = gn_dd_data_bcc.gn_dd_data_bcc()
    drv = gn_dd_data_hcp.gn_dd_data_hcp()

    dispatcher = {'hcp': drv.write_hcp_straight_data,
                  'hcporowan': drv.write_hcp_orawan_data,
                  'hcpprec': drv.inplane_hcp_beta1_prec,
                  'hcpten': drv.write_hcp_tensile_data,
                  'bccscrew': bcc.write_straight_screw_data}
    dispatcher[options.mtype.lower()]()
