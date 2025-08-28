from rxnrecer.cli import predict as rxnrecer_cli

def predict(data, mode, batch_size):
    res = rxnrecer_cli.step_by_step_prediction(
        input_data=data, 
        mode=mode,
        batch_size=batch_size
    )
    
    return res


if __name__ == "__main__":
    predict()
