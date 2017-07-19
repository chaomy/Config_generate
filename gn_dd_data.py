#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2017-07-19 11:32:14


from optparse import OptionParser
import gn_dd_data_hcp

if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype",
                      action="store",
                      type="string",
                      dest="mtype")
    (options, args) = parser.parse_args()

    drv = gn_dd_data_hcp.gn_dd_data_hcp()
    dispatcher = {'hcpdata': drv.write_hcp_straight_data,
                  'hcporowan': drv.write_hcp_orawan_data,
                  'hcpprec': drv.inplane_hcp_beta1_prec}
    dispatcher[options.mtype.lower()]()
