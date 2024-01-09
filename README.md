# Metodi Numerici Avanzati per Equazioni alle Derivate Parziali

Codici per l'esame del corso **Metodi Numerici Avanzati per Equazioni alle Derivate Parziali**.

## Codice

Implementazione dell'elemento di Crouzeix-Raviart per il problema di Stokes con dati al bordo di Dirichlet.

## Contenuto

- `/src/*` Codice sorgente per il FEM.
	- `/src/solver.m` Risolutore del prolema di Stokes sulla mesh.
    - `/src/basis.m` Funzioni di base per la velocità.
    - `/src/gradbasis.m` Gradienti delle funzioni di base per la velocità.
    - `/src/pressurebasis.m` Funzioni di base per la pressione.
    - `/src/loading.m` Sorgente.
    - `/src/exact.m` Soluzione esatta per la sorgente assegnata.
    - `/src/quadrature.m` Nodi di quadratura per l'elemento di riferimento.
    - `/src/errorsPlot.m` Grafici per l'andamento degli errori.
    - `/src/quadmeshes.mat` Varie mesh.