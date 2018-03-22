#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-10 23:05:38

import numpy as np


class gn_dd_load(object):

    def cal_stress(self, scale=7e6):
        burgs = self.burgs
        bkey = self.bkey
        burg = np.mat(burgs[bkey]['b'])
        norm = np.mat(burgs[bkey]['norms'][self.pln])
        print("burg", burg)
        print("norm", norm)
        # burgs = burgs.assign(e=e.values)
        # sigma = b /O n + n /O b
        sigma = np.outer(burg, norm) + np.outer(norm, burg)
        sigma = sigma * scale
        print(np.cross(burg, burg * sigma))

        # [sigma11, sigma22, sigma33, sigma23, sigma31, sigma12]
        applied = [sigma[0, 0], sigma[1, 1], sigma[2, 2],
                   sigma[1, 2], sigma[2, 0], sigma[0, 1]]

        print(applied)
        return applied

if __name__ == '__main__':
    drv = gn_dd_load()
    drv.cal_stress()
