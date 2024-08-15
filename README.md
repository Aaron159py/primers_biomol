# primers_biomol
Este código aun no está terminado
Me falta aplicarlo a secuencias y el método no_hairpin aun no esta terminado porque es un poco ambiguo

Para la selección de primers al momento de amplificar una secuencia deseada, se deben seguir ciertos criterios:
- Largo mínimo 17 bp
- No se deben formar primer-dimers; No más de 2 nucleótidos complementados en los últimos 4 nucleótidos del 3' end de los dos primers (forward y reverse)
- Al menos 50% GC content en ambos primers
- Extremo 3' G o C en cada primer para estabilidad
- T anneling 50-65 °C en ambos primers
- Diferencia máxima entre los primers de 5°C
- Ambos nucleotidos no deben formar hairpin; no mas de 2 nucleótidos haciendo hairpin (no entiendo muy bien esto)
