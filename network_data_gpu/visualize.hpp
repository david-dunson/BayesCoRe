
void plotImage(afarray A)
{
    af::Window window((int)A.dims(0), (int)A.dims(1), "2D plot example title");
    do
    {
        window.image(A);
        // window.plot(vertex_idx, S);
    } while (!window.close());
}

void plotHist(afarray A, int bins = 30)
{
    af::Window window(600, 600, "2D plot example title");

    afarray hist_out = histogram(A, bins);
    do
    {
        window.hist(hist_out, 0, 255);
        // window.plot(vertex_idx, S);
    } while (!window.close());
}