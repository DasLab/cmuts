/* preset.c -- the alignment presets and what each one expects.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "preset.h"

#include <stddef.h>

const cli_choice PRESET_CHOICES[] = {
    { "sr",       PRESET_SR,       "Short Read"                  },
    { "map-ont",  PRESET_MAP_ONT,  "Nanopore"                    },
    { "map-hifi", PRESET_MAP_HIFI, "PacBio HiFi"                 },
    { "map-pb",   PRESET_MAP_PB,   "PacBio CLR"                  },
    { "map-iclr", PRESET_MAP_ICLR, "Illumina Complete Long Read" },
    { "lr:hq",    PRESET_LR_HQ,    "Accurate Long Read"          },
    { NULL,       0,               NULL                          },
};

const char *preset_name(preset_id id)
{
    for (const cli_choice *choice = PRESET_CHOICES; choice->name; choice++) {
        if (choice->value == (int)id) {
            return choice->name;
        }
    }

    return NULL;
}

bool preset_accepts_pairs(preset_id id)
{
    switch (id) {
        case PRESET_SR:
            return true;
        case PRESET_MAP_ONT:
        case PRESET_MAP_HIFI:
        case PRESET_MAP_PB:
        case PRESET_MAP_ICLR:
        case PRESET_LR_HQ:
        case PRESET_UNSET:
            return false;
    }

    return false;
}
