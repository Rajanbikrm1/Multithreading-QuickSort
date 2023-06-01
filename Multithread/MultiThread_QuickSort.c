/**
 * @file lab2.c file contains all required code for Lab 2
 * @author Rajan Bikram Sah
 *
 * EECS 3540 – Operating Systems & Systems Programming
 * Spring 2023 – Programming Lab #2 – Multithreading Applications
 *
 * This project implements multi-thread Quicksort for sorting huge lists of integers.
 */

#define _GNU_SOURCE
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <stdbool.h>
#include <sys/time.h>
#include <sys/timeb.h>
#include <sys/resource.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <unistd.h>

// Global variables
bool finished[] = {false, false, false, false, false}; // An array of boolean values used to track whether each thread has finished or not.
int *numArray;                                         // A global variable to store the array that needs to be sorted.
int pivot[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};        // An array of integers used to store the pivot values.
int pivotArr[9];                                       // An array of integers used to store the pivot values for partitioning.

struct timeval start_load, end_load, start_sort, end_sort;
// Struct that stores the arguments for each thread.
struct thread_args
{
    int *arr;
    int low;
    int high;
    int early;
    int threshold;
    char alternate;
    char median;
};

/**
 * @brief A function to swap two integer values.
 *
 * @param a Pointer to the first integer value.
 * @param b Pointer to the second integer value.
 */
void swap(int *a, int *b)
{
    int temp = *a;
    *a = *b;
    *b = temp;
}

/**
 * @brief A function find the minimum value
 *
 * @param a the first integer value.
 * @param b the second integer value.
 */
int min(int a, int b)
{
    return (a < b) ? a : b;
}

// defining the load_array function 
void load_array(int *array, int size, int seed);
// defining the insertion_sort function 
void insertion_sort(int arr[], int low, int high);
// defining the shell_sort function
void shell_sort(int arr[], int low, int high);
// defining the quicksort function
void quicksort_single_thread(int *arr, int low, int high, int threshold, char alternate, char median);
// defining the quicksort function for multithreading
void *quicksort_multithread(void *args);
// defining the recursive_partition function for multithreading
void recursive_partition(int *array, int low, int high, int pieces, int *pivotArr, int *current_pieces);
// defining the findLargest function as helper function
int findLargest(int *array, int size);
// defining the thresholdCount function as helper function
int thresholdCount(int *numArray, int low, int high);
// defining the partition function as helper function for recursive_partition
int partition(int *numArray, int low, int high);
// defining the find_early_partition for early partition
int find_early_partition(int *array, int low, int high);
// defining the second_of_ten_partition for early implementation
int second_of_ten_partition(int *arr, int size);
// definng the function isSorted to check if the array is sorted or not
bool isSorted(int *numArray, int size);

/**
 * This is the main function of a program that performs quicksort on an array using
 * multiple threads. It receives several options through command-line arguments to
 * control its behavior, such as the array size, the seed for random number generation,
 * and the number of threads to use. It also measures the execution time of each stage
 * of the algorithm and prints statistics to the console at the end.
 *
 * @param argc The number of command-line arguments
 * @param argv An array of strings representing the command-line arguments
 *
 * @return 0 if the program completes successfully, or an error code otherwise
 */
int main(int argc, char *argv[])
{
    // Initialize default values for the program options
    int size = -1, threshold = 10, seed = -1, pieces = 10, max_threads = 4;
    char alternate = 'S', multithread = 'Y', median = 'N', early = 'N';

    // Declare variables for measuring the execution time of the program
    struct timeval start_time, end_time, start_sorting_time, end_sorting_time;
    struct rusage usage;
    double loading_time, sorting_time, wall_clock_time, cpu_time;

    // Parse the command-line options and set the corresponding variables
    int opt;
    while ((opt = getopt(argc, argv, "n:a:s:r:m:p:t:m3:e:")) != -1)
    {
        switch (opt)
        {
        case 'n':
            size = atoi(optarg);
            break;
        case 'a':
            alternate = optarg[0];
            break;
        case 's':
            threshold = atoi(optarg);
            break;
        case 'r':
            seed = atoi(optarg);
            break;
        case 'm':
            multithread = optarg[0];
            break;
        case 'p':
            pieces = atoi(optarg);
            break;
        case 't':
            max_threads = atoi(optarg);
            break;
        case 'm3':
            median = optarg[0];
            break;
        case 'e':
            early = optarg[0];
            break;
        default:
            printf("Usage: %s -n SIZE [-a ALTERNATE] [-s THRESHOLD] [-r SEED] [-m MULTITHREAD] [-p PIECES] [-t MAXTHREADS] [-m3 MEDIAN] [-e EARLY]\n", argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    // Check that the program options are valid
    if (size < 1 || size > 1000000000 || pieces < 1 || max_threads < 1 || max_threads > pieces)
    {
        printf("Invalid parameters\n");
        exit(EXIT_FAILURE);
    }

    if (!(alternate == 'S' || alternate == 's' || alternate == 'I' || alternate == 'i') ||
        !(multithread == 'Y' || multithread == 'y' || multithread == 'N' || multithread == 'n') ||
        !(median == 'Y' || median == 'y' || median == 'N' || median == 'n') ||
        !(early == 'Y' || early == 'y' || early == 'N' || early == 'n'))
    {
        printf("Invalid option values\n");
        exit(EXIT_FAILURE);
    }

    // Allocate memory for the array to be sorted
    numArray = (int *)malloc(size * sizeof(int));
    if (numArray == NULL)
    {
        fprintf(stderr, "Error allocating memory for array\n");
        exit(EXIT_FAILURE);
    }

    gettimeofday(&start_time, NULL);
    load_array(numArray, size, seed);       //load data from random.dat file into the numArray
    gettimeofday(&end_time, NULL);
    loading_time = (end_time.tv_sec - start_time.tv_sec) * 1000.0 + (end_time.tv_usec - start_time.tv_usec) / 1000.0;

    // printf("printing elements of array: \n");
    // for (int i = 0; i < size; i++)
    // {
    //     printf("%d ", numArray[i]);
    // }
    // printf("\n");

    // If early partition is enabled, create a thread for it
    pthread_t early_tid;
    if (early == 'Y' || early == 'y')
    {
        int early_partition = second_of_ten_partition(numArray, size);

        struct thread_args early_thread;
        pthread_attr_t early_attr;

        early_thread.early = 1;
        early_thread.low = 0;
        early_thread.high = early_partition;
        early_thread.median = median;

        pthread_attr_init(&early_attr);
        pthread_create(&early_tid, &early_attr, quicksort_multithread, (void *)&early_thread);
    }

    gettimeofday(&start_sorting_time, NULL);
    // Perform quicksort on the array using multiple threads if multithreading is enabled
    if (multithread == 'Y' || multithread == 'y')
    {

        // partition the array in 'pieces' pieces
        int pivotArr[pieces - 1];
        int current_pieces = 0;
        recursive_partition(numArray, 0, size - 1, pieces, pivotArr, &current_pieces);

        int nThreads = min(max_threads, pieces);
        struct thread_args threads[nThreads];
        pthread_t tid[nThreads];
        pthread_attr_t attr[nThreads];
        int thread_created = 0;

        // Create and run the initial threads
        for (int i = 0; i < nThreads; i++)
        {
            if (i == 0)
            {
                threads[i].low = 0;
                threads[i].high = pivotArr[i];
            }
            else if (i == pieces - 1)
            {
                threads[i].low = pivotArr[i - 1] + 1;
                threads[i].high = size - 1;
            }
            else
            {
                threads[i].low = pivotArr[i - 1] + 1;
                threads[i].high = pivotArr[i];
            }

            threads[i].threshold = threshold; // Pass the threshold value to the thread_args struct
            threads[i].arr = numArray;        // Set the arr field of the thread_args struct

            pthread_attr_init(&attr[i]);
            threads[i].median = median;
            pthread_create(&tid[i], &attr[i], quicksort_multithread, (void *)&threads[i]);
            thread_created++;
        }
        // Wait for the threads to finish
        for (int i = 0; i < thread_created; i++)
        {
            pthread_join(tid[i], NULL);
        }

        // The loop for managing additional partitions
        int remaining_partitions = pieces - nThreads;
        while (remaining_partitions > 0)
        {
            for (int i = 0; i < nThreads; i++)
            {
                if (pthread_tryjoin_np(tid[i], NULL) == 0)
                {
                    remaining_partitions--;
                    if (remaining_partitions > 0)
                    {
                        int part = thread_created % pieces;

                        if (part == 0)
                        {
                            threads[i].low = 0;
                            threads[i].high = pivotArr[0];
                        }
                        else if (part == pieces - 1)
                        {
                            threads[i].low = pivotArr[pieces - 2] + 1;
                            threads[i].high = size - 1;
                        }
                        else
                        {
                            threads[i].low = pivotArr[part - 1] + 1;
                            threads[i].high = pivotArr[part];
                        }

                        threads[i].threshold = threshold;                     // Pass the threshold value to the thread_args struct
                        threads[i].median = median;                           // Pass the median value to the thread_args struct
                        pthread_attr_init(&attr[i]);                          // initialize the new thread
                        pthread_create(&tid[i], &attr[i], quicksort_multithread, (void *)&threads[i]);
                        thread_created++; // Increment the total number of created threads

                        if (early == 'Y' || early == 'y')
                        {
                            // Find the early partition index
                            int index = find_early_partition(pivotArr, threads[i].low, threads[i].high);
                            if (index != -1)
                            {
                                pivotArr[index] = 0;
                            }
                        }
                    }
                }
            }
            // if none of the thread have finished,
            usleep(50000); // sleep 50 ms before checking again
        }
    }

    // If multithreading is disabled, perform single-threaded quicksort
    else
    {
        quicksort_single_thread(numArray, 0, size - 1, threshold, alternate, median);
    }

    // Wait for the early thread to finish if it was created
    if (early == 'Y' || early == 'y')
    {
        printf("Waiting for the early thread to finish...\n");
        pthread_join(early_tid, NULL);
    }

    // Display the summary statistics
    gettimeofday(&end_sorting_time, NULL);
    getrusage(RUSAGE_SELF, &usage);
    sorting_time = (end_sorting_time.tv_sec - start_sorting_time.tv_sec) * 1000.0 + (end_sorting_time.tv_usec - start_sorting_time.tv_usec) / 1000.0;
    wall_clock_time = (end_sorting_time.tv_sec - start_time.tv_sec) * 1000.0 + (end_sorting_time.tv_usec - start_time.tv_usec) / 1000.0;
    cpu_time = (usage.ru_utime.tv_sec + usage.ru_stime.tv_sec) * 1000.0 + (usage.ru_utime.tv_usec + usage.ru_stime.tv_usec) / 1000.0;

    printf("Load:   %f    ", loading_time / 1000);
    printf("Sort (Wall/CPU):   %f ", sorting_time / 1000);
    printf("  /   %f    ", cpu_time / 1000);
    printf("Total:   %f   ", wall_clock_time / 1000);
    if (isSorted(numArray, size))
    {
        printf("Sorted\n");
    }
    else
    {
        printf("Not sorted\n");
    }
    return 0;
}

/**
 * Loads a portion of integers from a binary file into an array, based on a given starting position and size.
 *
 * @param array Pointer to the array that will hold the loaded integers.
 * @param size The number of integers to load into the array.
 * @param seed The starting position in the binary file, or a negative value to generate a random seed.
 *
 * @throws std::runtime_error If the file cannot be opened.
 */
void load_array(int *array, int size, int seed)
{
    struct timeval start_time, end_time;
    gettimeofday(&start_time, NULL);

    FILE *file = fopen("random.dat", "rb"); // Open the binary file for reading
    // Check if file opened successfully
    if (file == NULL)
    {
        perror("Error opening file"); // Print an error message and exit the program if the file cannot be opened
        exit(EXIT_FAILURE);
    }

    if (seed < 0)
    {
        seed = rand() % 1000000000;
    }

    int start = seed;                           // Calculate the starting position in the file based on the seed value
    fseek(file, start * sizeof(int), SEEK_SET);     // Move the file pointer to the starting position
    int remaining = 600000000 - start; // Calculate the number of integers remaining in the file after the starting position
    
    //  Read integers from the file and load them into the array
    if (size <= remaining)
    {
        fread(array, sizeof(int), size, file);
    }
    else
    {
        fread(array, sizeof(int), remaining, file);
        fseek(file, 0, SEEK_SET); // Move the file pointer back to the beginning of the file
        fread(array + remaining, sizeof(int), size - remaining, file);
    }
    // Close the file
    fclose(file);
    gettimeofday(&end_time, NULL);
}

/**
 * Sorts a portion of an array in place, using the insertion sort algorithm.
 *
 * @param arr The array to be sorted.
 * @param low The index of the first element to be sorted.
 * @param high The index of the last element to be sorted.
 */
void insertion_sort(int arr[], int low, int high)
{
    for (int i = low + 1; i <= high; i++)
    {
        int key = arr[i];
        int j = i - 1;
        while (j >= low && arr[j] > key)
        {
            arr[j + 1] = arr[j];
            j--;
        }
        arr[j + 1] = key;
    }
}

/**
 * Sorts a portion of an array in place, using the shell sort algorithm.
 *
 * @param arr The array to be sorted.
 * @param low The index of the first element to be sorted.
 * @param high The index of the last element to be sorted.
 */
void shell_sort(int arr[], int low, int high)
{
    int n = high - low + 1;
    int h = 1;
    // Determine initial value of h
    while (h < n / 2)
    {
        h = 2 * h + 1;
    }
    // Keep decreasing value of h until it reaches 1
    while (h >= 1)
    {
        for (int i = h + low; i <= high; i++)
        {
            // Perform insertion sort on each segment of size h
            for (int j = i; j >= h + low && arr[j] < arr[j - h]; j -= h)
            {
                swap(&arr[j], &arr[j - h]);
            }
        }
        h = h / 2;
    }
}

/**
 * Sorts a portion of an array in place, using the quicksort algorithm, in a single thread.
 *
 * @param arr The array to be sorted.
 * @param low The index of the first element to be sorted.
 * @param high The index of the last element to be sorted.
 * @param threshold The threshold size for switching to another sorting algorithm.
 * @param alternate The alternative sorting algorithm to use for small partitions ('I' for insertion sort, 'S' for shell sort).
 * @param median Whether to use the median-of-three method for pivot selection ('Y' for yes, 'N' for no).
 */
void quicksort_single_thread(int *arr, int low, int high, int threshold, char alternate, char median)
{
    int size = high - low + 1;

    if (size < 2)
    {
        return;
    }
    else if (size == 2)
    {
        if (arr[low] > arr[high])
        {
            swap(&arr[low], &arr[high]);
        }
        return;
    }
    else if (size > 2 && size <= threshold)
    {
        if (alternate == 'S' || alternate == 's')
        {
            shell_sort(arr, low, high);
        }
        else if (alternate == 'I' || alternate == 'i')
        {
            insertion_sort(arr, low, high);
        }
        return;
    }
    else
    {
        if (median == 'Y' || median == 'y')
        {
            printf("Using median-of-three method\n");
            int mid = low + (high - low) / 2;
            if (arr[low] > arr[mid])
            {
                swap(&arr[low], &arr[mid]);
            }
            if (arr[low] > arr[high])
            {
                swap(&arr[low], &arr[high]);
            }
            if (arr[mid] > arr[high])
            {
                swap(&arr[mid], &arr[high]);
            }
            swap(&arr[low], &arr[mid]);
        }

        int pivot = arr[low];
        int i = low + 1;
        int j = high;

        while (1)
        {
            while (i <= j && arr[i] <= pivot)
                i++;

            while (arr[j] > pivot)
                j--;

            if (i <= j)
            {
                swap(&arr[i], &arr[j]);
            }
            else
            {
                break;
            }
        }

        swap(&arr[low], &arr[j]);
        if (j - low < high - i)
        {
            quicksort_single_thread(arr, low, j, threshold, alternate, median);
            quicksort_single_thread(arr, i, high, threshold, alternate, median);
        }
        else
        {
            quicksort_single_thread(arr, i, high, threshold, alternate, median);
            quicksort_single_thread(arr, low, j, threshold, alternate, median);
        }
    }
}

/**
 * Sorts a portion of an array in place, using the quicksort algorithm, in multiple threads.
 *
 * @param args A pointer to a struct containing the arguments for the function.
 */
void *quicksort_multithread(void *args)
{
    struct thread_args *t_args = (struct thread_args *)args;
    int *arr = t_args->arr;
    int low = t_args->low;
    int high = t_args->high;
    int threshold = t_args->threshold;
    char alternate = t_args->alternate;
    char median = t_args->median;

    quicksort_single_thread(arr, low, high, threshold, alternate, median);
    pthread_exit(NULL);
}

/**
 * Recursively partition an array into pieces using the QuickSort algorithm.
 *
 * @param array The array to partition.
 * @param low The lowest index of the sub-array to partition.
 * @param high The highest index of the sub-array to partition.
 * @param pieces The desired number of partitions.
 * @param pivotArr An array to store the partition points.
 * @param current_pieces The number of partitions currently created.
 *
 * @return None.
 */
void recursive_partition(int *array, int low, int high, int pieces, int *pivotArr, int *current_pieces)
{
    if (*current_pieces == pieces)
    {
        return;
    }

    int pi = partition(array, low, high);

    if (pi - 1 - low >= high - pi)
    {
        int new_high = pi - 1;
        int new_low = pi + 1;
        if (*current_pieces < pieces)
        {
            pivotArr[*current_pieces] = new_high;
            (*current_pieces)++;
        }
        recursive_partition(array, low, new_high, pieces, pivotArr, current_pieces);
        recursive_partition(array, new_low, high, pieces, pivotArr, current_pieces);
    }
    else
    {
        int new_high = high;
        int new_low = pi + 1;
        if (*current_pieces < pieces)
        {
            pivotArr[*current_pieces] = new_high;
            (*current_pieces)++;
        }
        recursive_partition(array, low, pi - 1, pieces, pivotArr, current_pieces);
        recursive_partition(array, new_low, high, pieces, pivotArr, current_pieces);
    }
}

/**
 * Finds the index of the largest partition in an array of partition points.
 *
 * @param arr The array of partition points.
 * @param size The size of the array.
 *
 * @return The index of the largest partition.
 */
int findLargest(int *arr, int size)
{
    int maxIndex = 0;
    for (int i = 1; i < size; i++)
    {
        if (arr[i] > arr[maxIndex])
        {
            maxIndex = i;
        }
    }
    return maxIndex;
}

/**
 * Counts the number of elements in a sub-array.
 *
 * @param numArray The array to count elements from.
 * @param low The lowest index of the sub-array.
 * @param high The highest index of the sub-array.
 *
 * @return The number of elements in the sub-array.
 */
int thresholdCount(int *numArray, int low, int high)
{
    int count = 0;
    for (int i = low; i <= high; i++)
    { // goes through loop and increment te size
        count++;
    }
    return count; // return the size
}

/**
 * Partitions an array using the QuickSort algorithm.
 *
 * @param numArray The array to partition.
 * @param low The lowest index of the sub-array to partition.
 * @param high The highest index of the sub-array to partition.
 *
 * @return The index of the partition point.
 */
int partition(int *numArray, int low, int high)
{
    int pivot = numArray[high]; // set the pivot at the last
    int i = (low - 1);          // set i before low
    for (int j = low; j <= high - 1; j++)
    { // goes through loop and arrange the order ascendingly
        if (numArray[j] <= pivot)
        { // arranging the order of the partition
            i++;
            int temp = numArray[i]; // interchanging the values
            numArray[i] = numArray[j];
            numArray[j] = temp;
        }
    }
    int temp2 = numArray[i + 1]; // interchanging the values with pivot
    numArray[i + 1] = numArray[high];
    numArray[high] = temp2;
    return (i + 1); // return the partiton point in the array
}

/**
 * Finds the index of an early partition point in an array of partition points.
 *
 * @param array The array of partition points.
 * @param low The index of the first partition point.
 * @param high The index of the last partition point.
 *
 * @return The index of the early partition point.
 */
int find_early_partition(int *array, int low, int high)
{
    for (int i = 0; i < 10; i++)
    {
        if (pivotArr[i] >= low && pivotArr[i] <= high)
        {
            return i;
        }
    }
    return -1;
}

/**
 * Partitions an array using the QuickSort algorithm, but with the pivot selected from the 2nd of 10 equally-spaced values.
 *
 * @param array The array to partition.
 * @param size The size of the array.
 *
 * @return The index of the partition point.
 */
int second_of_ten_partition(int *array, int size)
{
    int locations[11];
    int values[11];
    int step = size / 10;

    for (int i = 0; i < 11; i++)
    {
        locations[i] = i * step;
        values[i] = array[locations[i]];
    }

    int min_index = 0;
    // int second_min_index = 0;
    for (int i = 1; i < 11; i++)
    {
        if (values[i] < values[min_index])
        {
            min_index = i;
        }
    }
    values[min_index] = INT_MAX;
    int second_min_index = 0;
    for (int i = 1; i < 11; i++)
    {
        if (values[i] < values[second_min_index])
        {
            second_min_index = i;
        }
    }
    int x = locations[second_min_index];
    swap(&array[x], &array[0]);

    return partition(array, 0, size - 1);
}

/**
 * Checks if an array is sorted in ascending order.
 *
 * @param arr The array to be checked.
 * @param size The size of the array.
 *
 * @return 1 if the array is sorted, 0 otherwise.
 */
bool isSorted(int *numArray, int size)
{
    if (size == 1 || size == 0)
        return true; // if only one element then return true
    for (int i = 1; i < size; i++)
    { // goes in the loop and checks every corresponding values
        if (numArray[i - 1] > numArray[i])
            return false; // if finds one smaller returns false
    }
    return true; // otherwise returns true
}
