# OVL474---Tratamento_de_dados_astronomicos

### Photometric optimization
We added photometric error using LsstErrorModel from [photerr package (Crenshaw, J. et al., 2024)](https://github.com/jfcrenshaw/photerr) optimizing `nYrObs` parameter for each grizy band based on HSC Y3 data.

HSC Y3 query to obtain unique skymap_id and split our query into 10 chunks due to data volume and limit query time:

```sql
SELECT DISTINCT t.skymap_id
FROM pdr3_wide.patch_qa as t
ORDER BY t.skymap_id
```
HSC Y3 queries by ranges of skymap_id:
```sql
--- Version with selection by skymap values
WITH patches AS (
  SELECT skymap_id
  FROM pdr3_wide.patch_qa
  WHERE skymap_id BETWEEN 82780002 AND 91140304 --- Change skymap_id for each query
)
SELECT
  b.object_id, b.skymap_id, b.a_i,
  b.g_cmodel_mag, b.g_cmodel_magerr,
  b.r_cmodel_mag, b.r_cmodel_magerr,
  b.i_cmodel_mag, b.i_cmodel_magerr,
  b.z_cmodel_mag, b.z_cmodel_magerr,
  b.y_cmodel_mag, b.y_cmodel_magerr,
  b.a_g, b.a_r, b.a_z, b.a_y,
  b.ra, b.dec
FROM
  pdr3_wide.summary AS b
JOIN patches USING (skymap_id)
WHERE
  -- Basic flag cuts - Table 2 of Li et al. (2022) arXiv:2107.00136
  b.i_is_clean_centerpixels = True
  AND b.i_is_clean_allpixels = True
  AND b.i_sdsscentroid_flag = False
  AND b.i_extendedness_value != 0
  -- Galaxy property cuts
  AND b.i_cmodel_mag - b.a_i <= 24.5
  AND b.i_apertureflux_10_mag <= 25.5
```
Note that we are using similar cut to those used in [Li et al., 2022](https://arxiv.org/abs/2107.00136) but without shape-related cuts.
