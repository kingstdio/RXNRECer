def getUnirep(sequences, batch_size=10):
    """
    This function computes the unirep (h_avg) for each sequence in the input dataframe.
    The sequences are processed in batches to improve efficiency.
    
    Args:
        sequences (pd.DataFrame): DataFrame with columns 'id' and 'seq'.
        batch_size (int): The batch size for processing sequences.

    Returns:
        pd.DataFrame: DataFrame with 'id' and corresponding 'unirep' (h_avg).
    """
    # Initialize an empty list to store the batch results
    result_batches = []

    # Process sequences in batches
    for start in tqdm(range(0, len(sequences), batch_size)):
        end = min(start + batch_size, len(sequences))
        batch = sequences.iloc[start:end]
        
        # Zip id and seq together for the current batch
        batch_ids = batch['id'].tolist()
        batch_seqs = batch['seq'].tolist()

        # Call get_reps for the current batch
        h_avg, _, _ = get_reps(batch_seqs)
        
        # Combine the results (id, h_avg) for the current batch
        batch_result = pd.DataFrame({
            'id': batch_ids,
            'unirep': [h_avg[i] for i in range(len(h_avg))]
        })

        # Append the current batch result to the list
        result_batches.append(batch_result)

    # Concatenate all batch results into a single DataFrame
    final_result = pd.concat(result_batches, ignore_index=True)
    
    return final_result