# Metodi Numerici Avanzati per Equazioni alle Derivate Parziali

Codici per l'esame del corso **Metodi Numerici Avanzati per Equazioni alle Derivate Parziali**.

Implementazione dell'Elemento di Crouzeix-Raviart per risolvere il problema di Stokes con dati al bordo di Dirichlet omogenei.

## Contenuto

- `/src/*`
    - `/src/solver.m` Implementazione del risolutore per il problema di Stokes sulla mesh.
    - `/src/basis.m` Funzioni di base per la velocità.
    - `/src/gradbasis.m` Calcolo dei gradienti delle funzioni di base per la velocità.
    - `/src/pressurebasis.m` Funzioni di base per la pressione.
    - `/src/loading.m` Definizione della sorgente.
    - `/src/exact.m` Implementazione della soluzione esatta per la sorgente specificata.
    - `/src/quadrature.m` Nodi di quadratura per l'elemento di riferimento.
    - `/src/errorTrend.m` Generazione di grafici per visualizzare l'andamento degli errori.
    - `/src/reload.m` Wrapper per la funzione errorTrend.
    - `/src/quadmeshes.mat` Mesh predefinite.