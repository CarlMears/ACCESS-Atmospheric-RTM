// Public API for the ACCESS atmospheric RTM
//
// NOTE: be sure this is kept in sync with access_rtm.f90, there is no automatic
// checking between the two!
#pragma once

#include <stdint.h>

struct Era5DailyData {
    int32_t num_time, num_lats, num_lons, num_levels;

    int32_t *levels;
    float *lats;
    float *lons;
    int32_t *time;

    float *temperature;
    float *relative_humidity;
    float *height;
    float *liquid_content;
    float *surface_pressure;
    float *surface_temperature;
    float *surface_relative_humidity;
    float *surface_height;
    float *columnar_water_vapor;
    float *columnar_cloud_liquid;
};

struct RtmDailyData {
    int32_t num_hour, num_lat, num_lon, num_freq;

    int32_t *hour;
    float *lat;
    float *lon;
    float *freq;
    float *eia;

    float *col_vapor;
    float *col_water;

    float *tran;
    float *tb_up;
    float *tb_down;
};

int32_t compute_rtm(const struct Era5DailyData *atmo_data, struct RtmDailyData *rtm_data);