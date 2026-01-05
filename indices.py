veg_index_dict = {
    'NDVI': {
        'formula': lambda nir, red: (nir - red) / (nir + red),
        'bands': {'nir': 800, 'red': 660},
        'doi': '10.1016/S0034-4257(99)00067-X'
    },
    'GNDVI': {
        'formula': lambda nir, green, eps=1e-6: (nir - green) / (nir + green + eps),
        'bands': {'nir': 800, 'green': 550},
        'doi': '10.1016/S0034-4257(00)00113-9'
    },
    'NDVI705': {
        'formula': lambda r750, r705, eps=1e-6: (r750 - r705) / (r750 + r705 + eps),
        'bands': {'r750': 750, 'r705': 705},
        'doi': '10.1016/S0034-4257(02)00010-X'
    },
    'EVI': {
        'formula': lambda nir, red, blue, eps=1e-6: 2.5 * (nir - red) / (nir + 6*red - 7.5*blue + 1 + eps),
        'bands': {'nir': 800, 'red': 660, 'blue': 470},
        'doi': '10.1016/S0034-4257(02)00096-2'
    },
    'OSAVI': {
        'formula': lambda nir, red, eps=1e-6: (nir - red) / (nir + red + 0.16 + eps),
        'bands': {'nir': 800, 'red': 660},
        'doi': '10.1016/0034-4257(95)00186-7'
    },
    'MCARI': {
        'formula': lambda nir, red, green, eps=1e-6: ((nir - red) - 0.2 * (nir - green)) * (nir / (red + eps)),
        'bands': {'nir': 700, 'red': 670, 'green': 550},
        'doi': '10.1016/S0034-4257(00)00113-9'
    },
    'TCARI': {
        'formula': lambda nir, red, green, eps=1e-6: 3 * ((nir - red) - 0.2 * (nir - green) * (nir / (red + eps))),
        'bands': {'nir': 700, 'red': 670, 'green': 550},
        'doi': '10.1016/S0034-4257(02)00018-4'
    },
    'CIgreen': {
        'formula': lambda nir, green, eps=1e-6: nir / (green + eps) - 1,
        'bands': {'nir': 800, 'green': 550},
        'doi': '10.1029/2002GL016450'
    },
    'CIred-edge': {
        'formula': lambda nir, red_edge, eps=1e-6: nir / (red_edge + eps) - 1,
        'bands': {'nir': 750, 'red_edge': 710},
        'doi': '10.1029/2005GL022688'
    },
    'PRI': {
        'formula': lambda r531, r570, eps=1e-6: (r531 - r570) / (r531 + r570 + eps),
        'bands': {'r531': 531, 'r570': 570},
        'doi': '10.1080/014311697217387'
    },
    'PSRI': {
        'formula': lambda r678, r500, r750, eps=1e-6: (r678 - r500) / (r750 + eps),
        'bands': {'r678': 678, 'r500': 500, 'r750': 750},
        'doi': '10.1034/j.1399-3054.1999.106119.x'
    },
    'MTCI': {
        'formula': lambda nir, red_edge, red, eps=1e-6: (nir - red_edge) / (red_edge - red + eps),
        'bands': {'nir': 754, 'red_edge': 709, 'red': 681},
        'doi': '10.1080/0143116042000274015'
    },
    'NPCI': {
        'formula': lambda r680, r430, eps=1e-6: (r680 - r430) / (r680 + r430 + eps),
        'bands': {'r680': 680, 'r430': 430},
        'doi': '10.1016/0034-4257(94)90136-8'
    },
    'WBI': {
        'formula': lambda r900, r970, eps=1e-6: r900 / (r970 + eps),
        'bands': {'r900': 900, 'r970': 970},
        'doi': '10.1016/S0034-4257(00)00147-4'
    }
}
